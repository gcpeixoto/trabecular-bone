import sys
import os
from subprocess import call

import xml.etree.ElementTree as etree
import xml.dom.minidom as minidom
# from lxml import etree

import numpy as np
from collections import OrderedDict

from read_xplt import ELEM_TYPE

import pdb


def el_type_number(feb_elem_type):

    if feb_elem_type == ELEM_TYPE['ELEM_HEX']:
        vtk_elem_type = 12
    elif feb_elem_type == ELEM_TYPE['ELEM_PENTA']:
        vtk_elem_type = 12
    elif feb_elem_type == ELEM_TYPE['ELEM_TET']:
        vtk_elem_type = 10
    elif feb_elem_type == ELEM_TYPE['ELEM_QUAD']:
        vtk_elem_type = 9
    elif feb_elem_type == ELEM_TYPE['ELEM_TRI']:
        vtk_elem_type = 5

    return vtk_elem_type


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """

    rough_string = etree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)

    return reparsed.toprettyxml(indent="\t")


def load_geom(workdir, ndfiles, elfiles):

    nodes = []
    elems = []

    # load data
    for ndf in ndfiles:
        nodes.append(np.loadtxt(workdir + ndf + '.dat'))
    for elf in elfiles:
        elems.append(np.loadtxt(workdir + elf + '.dat', dtype=int)[:, 1:])

    # append nodes
    if len(ndfiles) < 2:
        node_coords = nodes[0]
    else:
        for i in range(len(ndfiles) - 1):
            node_coords = np.vstack((nodes[i], nodes[i + 1]))

    dom_n_elems = map(len, elems)

    return node_coords, elems, dom_n_elems


def load_data(workdir, nstate, n_subdomains):

    # load miscelleneous data
    dom_elem_types = np.loadtxt(workdir + 'element_types_%d.dat' %
                                nstate, dtype=int, ndmin=1)

    item_formats = np.loadtxt(workdir + 'item_format_%d.dat' % nstate,
                              dtype=int, ndmin=1)

    item_names = np.genfromtxt(workdir + 'item_names_%d.dat' % nstate,
                               dtype='str', delimiter='\n')

    # load node data
    item_def_doms = []
    with open(workdir + 'item_def_doms_%d.dat' % nstate, 'r') as f:
        for line in f:
            item_def_doms.append(map(int, line[:-2].split(' ')))

    # assumption: node data is define on all the nodes
    n_node_data = item_def_doms.count([-1])
    node_data = []
    for nf in item_names[:n_node_data]:
        filename = workdir + '%s_%d.dat' % (nf, nstate)
        node_data.append(np.loadtxt(filename))

    # load element data
    elem_files = []
    elem_data = []
    data_dims = [0] * (len(item_names) - n_node_data)
    k = 0
    for eldn in item_names[n_node_data:]:
        elf = []
        eld = []
        for i in range(n_subdomains):
            filename = workdir + '%s_%d_%d.dat' % (eldn, nstate, i)
            if os.path.isfile(filename):
                elf.append(filename)
                eld.append(np.loadtxt(filename))

                if eld[i].ndim == 1:
                    data_dims[k] = 1
                    eld[i] = np.array([eld[i]]).T
                else:
                    data_dims[k] = len(eld[i][0])

            else:
                elf.append('')
                eld.append([])

        elem_files.append(elf)
        elem_data.append(eld)
        k += 1

    return node_data, elem_data,\
        dom_elem_types, item_formats, item_names, item_def_doms, data_dims


def write_vtk(workdir, nodes, elems, dom_n_elems, node_data,
              elem_data, dom_elem_types, item_formats, item_names,
              item_def_doms, data_dims, vtkfile):

    output_file = workdir + vtkfile

    # material type (subdomain number)
    mat_types = []
    for i in range(len(dom_n_elems)):
        mat_types += [i] * dom_n_elems[i]
    mat_types = np.array(mat_types, dtype=int)

    # OUTPUT
    print "Outputting to %s..." % output_file

    # header
    root = etree.Element("VTKFile", type="UnstructuredGrid")

    UnstructuredGrid_xml = etree.SubElement(root, "UnstructuredGrid")

    Piece_xml = etree.SubElement(UnstructuredGrid_xml, "Piece",
                                 NumberOfPoints=str(len(nodes)),
                                 NumberOfCells=str(sum(dom_n_elems)))

    PointData_xml = etree.SubElement(Piece_xml, "PointData")

    for i in range(len(node_data)):
        if (type(node_data[i][0]) is np.ndarray):
            data_dim = len(node_data[i][0])
        elif (type(node_data[i][0]) is np.float64):
            data_dim = 1
        else:
            assert(0), "Wrong data type"

        DataArray_xml = etree.SubElement(PointData_xml, "DataArray",
                                         type="Float32",
                                         NumberOfComponents=str(data_dim),
                                         Name=item_names[i])
        DataArray_xml.text = '\n' + ' '.join([repr(a) for a in
                                              node_data[i].flatten()]) + '\n'

    # for each node parameter
    for i in range(len(elem_data)):
        if item_formats[i + len(node_data)] == 0:
            def_doms = item_def_doms[i + len(node_data)]

            # skip for undefined parameter
            if def_doms == [-2]:
                continue

            nodes_set = set([])
            # fill zeros to the empty data
            node_data_all = np.zeros([len(nodes), data_dims[i]])

            # gather all the 'defined' nodal values from the
            # subdomains
            nodes_set = []
            for j in def_doms:
                nodes_set += list(elems[j].flatten())

                # remove duplicates while preserving order
                nodes_set = list(OrderedDict.fromkeys(nodes_set))
                for k in range(len(nodes_set)):
                    node_data_all[nodes_set[k]] = elem_data[i][j][k]

            DataArray_xml \
                = etree.SubElement(PointData_xml, "DataArray",
                                   type="Float32",
                                   NumberOfComponents=str(data_dims[i]),
                                   Name=item_names[i + len(node_data)])
            DataArray_xml.text \
                = '\n' + ' '.join([repr(a) for a in
                                   node_data_all.flatten()]) + '\n'

    CellData_xml = etree.SubElement(Piece_xml, "CellData")

    # material types
    DataArray_xml = etree.SubElement(CellData_xml, "DataArray",
                                     type="Int32",
                                     NumberOfComponents="1",
                                     Name="mat_types")
    DataArray_xml.text = '\n' + ' '.join([`a` for a in mat_types]) + '\n'

    # for each element parameter
    for i in range(len(elem_data)):
        if item_formats[i + len(node_data)] == 1:
            def_doms = item_def_doms[i + len(node_data)]
            # skip for undefined parameter
            if def_doms == [-2]:
                continue

            # gather all the 'defined' element values from the
            # subdomains
            # fill zeros to the empty data
            for i_ed in range(len(elem_data[i])):

                if i_ed not in def_doms:
                    elem_data[i][i_ed] = np.zeros([dom_n_elems[i_ed],
                                                   data_dims[i]])

                data_type = str(elem_data[i][i_ed].dtype)
                if 'int' in data_type:
                    data_type = "Int64"
                elif 'float' in data_type:
                    data_type = "Float64"
                else:
                    data_type = "Float64"

            eld_txt = ''
            for i_ed in range(len(elem_data[i])):
                eld_txt += '\n' + ' '.join([repr(a) for a in
                                            elem_data[i][i_ed].flatten()]) + '\n'

            DataArray_xml \
                = etree.SubElement(CellData_xml, "DataArray",
                                   type=data_type,
                                   NumberOfComponents=str(data_dims[i]),
                                   Name=item_names[i + len(node_data)])
            DataArray_xml.text = eld_txt

    Points_xml = etree.SubElement(Piece_xml, "Points")
    DataArray_xml = etree.SubElement(Points_xml, "DataArray",
                                     NumberOfComponents="3",
                                     type="Float32")
    DataArray_xml.text = '\n' + ' '.join([`a` for a in nodes.flatten()]) + '\n'

    Cells_xml = etree.SubElement(Piece_xml, "Cells")
    DataArray_xml = etree.SubElement(Cells_xml, "DataArray",
                                     type="Int32",
                                     Name="connectivity")
    elem_conns = ''
    for i in range(len(elems)):
        elem_conns += '\n' + ' '.join([`a` for a in elems[i].flatten()]) + '\n'
    DataArray_xml.text = elem_conns

    DataArray_xml = etree.SubElement(Cells_xml, "DataArray",
                                     type="Int32",
                                     Name="offsets")
    elem_offs = []
    for i in range(len(elems)):
        elem_offs += [len(a) for a in elems[i]]
    elem_offs = np.cumsum(elem_offs)
    DataArray_xml.text = '\n' + ' '.join([`a` for a in elem_offs]) + '\n'

    DataArray_xml = etree.SubElement(Cells_xml, "DataArray",
                                     type="UInt8",
                                     Name="types")

    # assume that each subdomain only consists of the same element type
    dom_elem_types = map(el_type_number, dom_elem_types)
    elem_types = ''
    for i in range(len(dom_elem_types)):
        elem_types += '\n' + \
            ' '.join([str(dom_elem_types[i])] * dom_n_elems[i]) + '\n'
    DataArray_xml.text = elem_types

    # xml_str = etree.tostring(root, encoding='ISO-8859-1',
    # pretty_print="true")
    xml_str = prettify(root)

    # output
    f_out = open(output_file, 'w')
    f_out.write(xml_str)
    f_out.close()
