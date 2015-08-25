"""
read_xplt.py

Read FEBIO's xplt format and output the nodde/element data.

Ref: http://mrldata.sci.utah.edu/data/FEBioBinaryDatabaseSpecification.pdf

"""

import struct
from numpy import *
import sys
import pdb


TAGS = {
    'FEBIO': '0x00464542',
    'VERSION': '0x0004',
    'ROOT': '0x01000000',
    'HEADER': '0x01010000',
    'HDR_VERSION': '0x01010001',
    'HDR_NODES': '0x01010002',
    'HDR_MAX_FACET_NODES': '0x01010003',
    'HDR_COMPRESSION': '0x01010004',  # addedinversion0x03
    'DICTIONARY': '0x01020000',
    'DIC_ITEM': '0x01020001',
    'DIC_ITEM_TYPE': '0x01020002',
    'DIC_ITEM_FMT': '0x01020003',
    'DIC_ITEM_NAME': '0x01020004',
    'DIC_GLOBAL': '0x01021000',
    'DIC_MATERIAL': '0x01022000',
    'DIC_NODAL': '0x01023000',
    'DIC_DOMAIN': '0x01024000',
    'DIC_SURFACE': '0x01025000',
    'MATERIALS': '0x01030000',
    'MATERIAL': '0x01030001',
    'MAT_ID': '0x01030002',
    'MAT_NAME': '0x01030003',
    'GEOMETRY': '0x01040000',
    'NODE_SECTION': '0x01041000',
    'NODE_COORDS': '0x01041001',
    'DOMAIN_SECTION': '0x01042000',
    'DOMAIN': '0x01042100',
    'DOMAIN_HDR': '0x01042101',
    'DOM_ELEM_TYPE': '0x01042102',
    'DOM_MAT_ID': '0x01042103',
    'DOM_ELEMS': '0x01032104',
    'DOM_NAME': '0x01032105',  # addedinversion0x03
    'DOM_ELEM_LIST': '0x01042200',
    'ELEMENT': '0x01042201',
    'SURFACE_SECTION': '0x01043000',
    'SURFACE': '0x01043100',
    'SURFACE_HDR': '0x01043101',
    'SURFACE_ID': '0x01043102',
    'SURFACE_FACES': '0x01043103',
    'SURFACE_NAME': '0x01043104',  # addedinversion0x03
    'FACE_LIST': '0x01043200',
    'FACE': '0x01043201',
    'NODESET_SECTION': '0x01044000',  # addedinversion0x03
    'NODESET': '0x01044100',  # addedinversion0x03
    'NODESET_HDR': '0x01044101',  # addedinversion0x03
    'NODESET_ID': '0x01044102',  # addedinversion0x03
    'NODESET_NAME': '0x01044103',  # addedinversion0x03
    'NODESET_SIZE': '0x01044104',  # addedinversion0x03
    'NODESET_LIST': '0x01044200',  # addedinversion0x03
    'STATE': '0x02000000',
    'STATE_HEADER': '0x02010000',
    'STATE_HDR_ID': '0x02010001',
    'STATE_HDR_TIME': '0x02010002',
    'STATE_DATA': '0x02020000',
    'STATE_VARIABLE': '0x02020001',
    'STATE_VAR_ID': '0x02020002',
    'STATE_VAR_DATA': '0x02020003',
    'GLOBAL_DATA': '0x02020100',
    'MATERIAL_DATA': '0x02020200',
    'NODE_DATA': '0x02020300',
    'ELEMENT_DATA': '0x02020400',
    'FACE_DATA': '0x02020500'
}


ELEM_TYPE = {
    'ELEM_HEX': 0,
    'ELEM_PENTA': 1,
    'ELEM_TET': 2,
    'ELEM_QUAD': 3,
    'ELEM_TRI': 4,
    'ELEM_TRUSS': 5,
    'ELEM_HEX20': 6,
    'ELEM_TET10': 7,
    'ELEM_TET15': 8,
    'ELEM_HEX27': 8
}


def num_el_nodes(dom_elem_type, ELEM_TYPE):

    if dom_elem_type == ELEM_TYPE['ELEM_HEX']:
        ne = 8
    elif dom_elem_type == ELEM_TYPE['ELEM_PENTA']:
        ne = 6
    elif dom_elem_type == ELEM_TYPE['ELEM_TET']:
        ne = 4
    elif dom_elem_type == ELEM_TYPE['ELEM_QUAD']:
        ne = 4
    elif dom_elem_type == ELEM_TYPE['ELEM_TRI']:
        ne = 3
    elif dom_elem_type == ELEM_TYPE['ELEM_TRUSS']:
        ne = 2
    elif dom_elem_type == ELEM_TYPE['ELEM_HEX20']:
        ne = 20
    elif dom_elem_type == ELEM_TYPE['ELEM_TET10']:
        ne = 10
    elif dom_elem_type == ELEM_TYPE['ELEM_TET15']:
        ne = 15
    elif dom_elem_type == ELEM_TYPE['ELEM_HEX27']:
        ne = 27

    return ne


def search_block(f, TAGS, BLOCK_TAG, max_depth=5, cur_depth=0,
                 verbose=0, inv_TAGS=0, print_tag=0):
    '''Search some block in the current level.'''

    # record the initial cursor position
    if cur_depth == 0:
        ini_pos = f.tell()

    if cur_depth > max_depth:
        print 'Max iteration reached: Cannot find %s' % BLOCK_TAG
        return -1

    buf = f.read(4)

    if buf == '':
        print 'EOF: Cannot find %s' % BLOCK_TAG
        return -1
    else:
        cur_id = struct.unpack('I', buf)[0]
    a = struct.unpack('I', f.read(4))[0]  # size of the block

    if verbose == 1:
        cur_id_str = '0x' + '{0:08x}'.format(cur_id)
        # print 'cur_ID: ' + cur_id_str
        print 'cur_tag:', inv_TAGS[cur_id_str]
        print 'size:', a

    if(int(TAGS[BLOCK_TAG], base=16) == cur_id):
        if print_tag == 1:
            print('%s' % BLOCK_TAG)

        return a
    else:
        f.seek(a, 1)
        d = search_block(f, TAGS, BLOCK_TAG, cur_depth=cur_depth + 1,
                         verbose=verbose,
                         inv_TAGS=inv_TAGS,
                         print_tag=print_tag)
        if d == -1:
            # put the cursor position back
            if cur_depth == 0:
                f.seek(ini_pos, 0)
            return -1
        else:
            return d


def check_block(f, TAGS, BLOCK_TAG, filesize=-1):
    '''Check if the BLOCK TAG exists immediately after the file cursor.'''

    if filesize > 0:
        if f.tell() + 4 > filesize:
            print "EOF reached"
            return 0

    buf = struct.unpack('I', f.read(4))[0]
    f.seek(-4, 1)

    if(int(TAGS[BLOCK_TAG], base=16) == buf):
        return 1

    return 0


def seek_block(f, TAGS, BLOCK_TAG):
    if(int(TAGS[BLOCK_TAG], base=16) == struct.unpack('I', f.read(4))[0]):
        print('%s' % BLOCK_TAG)
    a = struct.unpack('I', f.read(4))  # size of the root section

    return a[0]


def read_xplt(workdir, filename, nstate, TAGS):

    # output files to be returned
    ndfiles = []
    elfiles = []

    if workdir[-1] is not '/':
        workdir += '/'

    f = open(workdir + filename, 'rb')

    inv_TAGS = {v: k for k, v in TAGS.items()}

    # get file size
    f.seek(0, 2)
    filesize = f.tell()
    f.seek(0, 0)

    if(int(TAGS['FEBIO'], base=16) == struct.unpack('I', f.read(4))[0]):
        print('Correct FEBio format')

    search_block(f, TAGS, 'ROOT')

    search_block(f, TAGS, 'HEADER')

    search_block(f, TAGS, 'HDR_VERSION')
    version = struct.unpack('I', f.read(4))[0]
    if(int(TAGS['VERSION'], base=16) == version):
        print('Current version is: %d' % version)
    else:
        print('Incorrect plot file version!')
        return -1

    a = search_block(f, TAGS, 'HDR_NODES')
    nNodes = struct.unpack('I', f.read(4))[0]
    print 'Number of nodes: %d' % (nNodes)

    search_block(f, TAGS, 'DICTIONARY')

    item_types = []
    item_formats = []  # 0: nodeal values, 1: elemental values
    item_names = []
    a = search_block(f, TAGS, 'DIC_NODAL')
    while check_block(f, TAGS, 'DIC_ITEM'):

        a = search_block(f, TAGS, 'DIC_ITEM')
        a = search_block(f, TAGS, 'DIC_ITEM_TYPE')
        item_types.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'DIC_ITEM_FMT')
        item_formats.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'DIC_ITEM_NAME')
        item_names.append(f.read(64).split('\x00')[0])

    a = search_block(f, TAGS, 'DIC_DOMAIN')
    while check_block(f, TAGS, 'DIC_ITEM'):

        a = search_block(f, TAGS, 'DIC_ITEM')
        a = search_block(f, TAGS, 'DIC_ITEM_TYPE')
        item_types.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'DIC_ITEM_FMT')
        item_formats.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'DIC_ITEM_NAME')
        item_names.append(f.read(64).split('\x00')[0])

    a = search_block(f, TAGS, 'MATERIALS')
    mat_names = []
    mat_ids = []
    while check_block(f, TAGS, 'MATERIAL'):
        a = search_block(f, TAGS, 'MATERIAL')

        a = search_block(f, TAGS, 'MAT_ID')
        mat_ids.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'MAT_NAME')
        mat_names.append(f.read(64).split('\x00')[0])

    a = search_block(f, TAGS, 'GEOMETRY')
    a = search_block(f, TAGS, 'NODE_SECTION')

    a = search_block(f, TAGS, 'NODE_COORDS')
    n_nodes = int(a / 3 / 4)
    node_coords = zeros([n_nodes, 3])
    for i in range(n_nodes):
        for j in range(0, 3):
            node_coords[i, j] = struct.unpack('f', f.read(4))[0]
    savetxt(workdir + 'nodes_%d.dat' % nstate, node_coords)
    ndfiles.append('nodes_%d' % nstate)

    a = search_block(f, TAGS, 'DOMAIN_SECTION')
    dom_elem_types = []
    dom_mat_ids = []
    dom_names = []
    dom_n_elems = []  # number of elements for each domain
    dom_elements = []  # elements for each domain
    # NOTE: index starts from 0 (in .feb file, index starts from 1)
    while check_block(f, TAGS, 'DOMAIN'):
        a = search_block(f, TAGS, 'DOMAIN')

        a = search_block(f, TAGS, 'DOMAIN_HDR')

        a = search_block(f, TAGS, 'DOM_ELEM_TYPE')
        dom_elem_type = int(struct.unpack('I', f.read(4))[0])
        dom_elem_types.append(dom_elem_type)

        a = search_block(f, TAGS, 'DOM_MAT_ID')
        dom_mat_ids.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'DOM_ELEMS')
        dom_n_elems.append(int(struct.unpack('I', f.read(4))[0]))

        # a = search_block(f, TAGS, 'DOM_NAME', verbose=1, inv_TAGS=inv_TAGS)
        # dom_names.append(int(struct.unpack('I', f.read(4))[0]))

        a = search_block(f, TAGS, 'DOM_ELEM_LIST')

        ne = num_el_nodes(dom_elem_type, ELEM_TYPE)
        elements = []
        while check_block(f, TAGS, 'ELEMENT'):
            a = search_block(f, TAGS, 'ELEMENT', print_tag=0)
            element = zeros(ne + 1, dtype=int)
            for j in range(ne + 1):
                element[j] = struct.unpack('I', f.read(4))[0]
            elements.append(element)

        dom_elements.append(elements)

    for i in range(len(dom_elements)):
        savetxt(workdir + 'elements_%d_%d.dat' % (nstate, i),
                dom_elements[i], fmt='%d')
        elfiles.append('elements_%d_%d' % (nstate, i))

    if search_block(f, TAGS, 'SURFACE_SECTION') > 0:

        surface_ids = []
        surface_faces = []  # number of faces
        surface_names = []
        faces = []
        face_ids = []
        while check_block(f, TAGS, 'SURFACE'):
            a = search_block(f, TAGS, 'SURFACE')

            a = search_block(f, TAGS, 'SURFACE_HDR')

            a = search_block(f, TAGS, 'SURFACE_ID')
            surface_ids.append(struct.unpack('I', f.read(4))[0])

            # number of facets
            a = search_block(f, TAGS, 'SURFACE_FACES')
            surface_faces.append(struct.unpack('I', f.read(4))[0])

            a = search_block(f, TAGS, 'SURFACE_NAME')
            # surface name length is specified just above
            surface_names.append(f.read(a).split('\x00')[0])

            if (check_block(f, TAGS, 'FACE_LIST') == 0):
                continue
            else:
                a = search_block(f, TAGS, 'FACE_LIST')

            while check_block(f, TAGS, 'FACE'):
                a = search_block(f, TAGS, 'FACE')
                cur_cur = f.tell()

                face = zeros(3, dtype=int)
                face_ids.append(struct.unpack('I', f.read(4))[0])

                # skip (probably specifing the surface element type here)
                f.seek(4, 1)
                # tri3 element
                for j in range(3):
                    face[j] = struct.unpack('I', f.read(4))[0]
                faces.append(face)
                # skip junk
                f.seek(cur_cur + a, 0)

    # skip nodeset section
    a = search_block(f, TAGS, 'NODESET_SECTION', print_tag=0)
    cur_cur = f.tell()
    f.seek(cur_cur + a, 0)

    # skip the first nstate states
    cur_state = 0
    while check_block(f, TAGS, 'STATE', filesize=filesize) &\
            (cur_state < nstate):
        a_state = search_block(f, TAGS, 'STATE', print_tag=0)

        cur_cur = f.tell()
        a = search_block(f, TAGS, 'STATE_HEADER', print_tag=0)
        # a = search_block(f, TAGS, 'STATE_HDR_ID', print_tag=0)
        a = search_block(f, TAGS, 'STATE_HDR_TIME', print_tag=0)
        time = struct.unpack('f', f.read(4))
        print 'This state is at %f time' % (time)
        f.seek(cur_cur + a_state, 0)

        # f.seek(a_state, 1)

        cur_state += 1

    if cur_state != nstate:
        print "State %d does not exist!" % nstate
        print "Current state is, ", cur_state
        return -1

    # now extract the information from the desired state
    a = search_block(f, TAGS, 'STATE')
    a = search_block(f, TAGS, 'STATE_HEADER')

    # a = search_block(f, TAGS, 'STATE_HDR_ID', verbose=1, inv_TAGS=inv_TAGS)
    # state_id = struct.unpack('f', f.read(4))
    # print 'The state ID: %d' % state_id

    a = search_block(f, TAGS, 'STATE_HDR_TIME')
    time = struct.unpack('f', f.read(4))
    print 'This state is at %f time' % (time)

    a = search_block(f, TAGS, 'STATE_DATA')

    # a = search_block(f, TAGS, 'MATERIAL_DATA', verbose=1, inv_TAGS=inv_TAGS)

    a = search_block(f, TAGS, 'NODE_DATA')
    n_node_data = 0
    item_def_doms = []
    while check_block(f, TAGS, 'STATE_VARIABLE'):
        n_node_data += 1

        a = search_block(f, TAGS, 'STATE_VARIABLE')
        a = search_block(f, TAGS, 'STATE_VAR_ID')
        var_id = struct.unpack('I', f.read(4))[0]
        print '\nvariable id:', var_id
        print 'variable_name:', item_names[var_id - 1]

        a = search_block(f, TAGS, 'STATE_VAR_DATA')

        a_end = f.tell() + a
        if item_types[var_id - 1] == 0:  # FLOAT
            data_dim = 1
        elif item_types[var_id - 1] == 1:  # VEC3F
            data_dim = 3
        # MAT3FS (6 elements due to symmetry)
        elif item_types[var_id - 1] == 2:
            data_dim = 6
        else:
            print 'unknwon data dimension!'
            return -1

        # assumption: node data is defined for all the ndoes
        def_doms = []
        while(f.tell() < a_end):
            dom_num = struct.unpack('I', f.read(4))[0]
            data_size = struct.unpack('I', f.read(4))[0]
            n_data = int(data_size / data_dim / 4.0)
            def_doms.append(dom_num - 1)

            print 'number of node data for domain %s = %d'\
                % (dom_num, n_data)

            if n_data > 0:
                elem_data = zeros([n_data, data_dim])
                for i in range(0, n_data):
                    for j in range(0, data_dim):
                        elem_data[i, j] = struct.unpack('f', f.read(4))[0]
                savetxt(workdir + '%s_%d.dat'
                        % (item_names[var_id - 1], nstate),
                        elem_data)

        item_def_doms.append(def_doms)

    a = search_block(f, TAGS, 'ELEMENT_DATA')
    while check_block(f, TAGS, 'STATE_VARIABLE'):

        a = search_block(f, TAGS, 'STATE_VARIABLE')

        a = search_block(f, TAGS, 'STATE_VAR_ID')
        var_id = struct.unpack('I', f.read(4))[0] + n_node_data
        print '\nvariable id:', var_id
        print 'variable_name:', item_names[var_id - 1]

        a = search_block(f, TAGS, 'STATE_VAR_DATA')
        a_end = f.tell() + a
        if item_types[var_id - 1] == 0:  # FLOAT
            data_dim = 1
        elif item_types[var_id - 1] == 1:  # VEC3F
            data_dim = 3
        # MAT3FS (6 elements due to symmetry)
        elif item_types[var_id - 1] == 2:
            data_dim = 6
        else:
            print 'unknwon data dimension!'
            return -1

        def_doms = []
        while(f.tell() < a_end):
            dom_num = struct.unpack('I', f.read(4))[0]
            data_size = struct.unpack('I', f.read(4))[0]
            n_data = int(data_size / data_dim / 4.0)
            def_doms.append(dom_num - 1)
            print 'number of element data for domain %s = %d'\
                % (dom_num, n_data)

            if n_data > 0:
                elem_data = zeros([n_data, data_dim])
                for i in range(0, n_data):
                    for j in range(0, data_dim):
                        elem_data[i, j] = struct.unpack('f', f.read(4))[0]
                savetxt(workdir + '%s_%d_%d.dat'
                        % (item_names[var_id - 1], nstate, dom_num - 1),
                        elem_data)

        item_def_doms.append(def_doms)

        if f.tell() >= filesize:
            break

    print '\n\ndata extraction done'

    # skip the last whatever states
    cur_state = nstate
    if f.tell() < filesize:
        while check_block(f, TAGS, 'STATE'):
            cur_state += 1
            a_state = search_block(f, TAGS, 'STATE', print_tag=0)

            cur_cur = f.tell()
            a = search_block(f, TAGS, 'STATE_HEADER', print_tag=0)
            a = search_block(f, TAGS, 'STATE_HDR_TIME', print_tag=0)
            time = struct.unpack('f', f.read(4))
            # print 'This state is at %f time' % (time)
            f.seek(cur_cur, 0)

            f.seek(a_state, 1)
            if f.tell() == filesize:
                break

    if f.tell() == filesize:
        print 'EOF reached.'

    # save element types and item formats for each subdomain
    savetxt(workdir + 'element_types_%d.dat' % (nstate),
            dom_elem_types, fmt='%d')
    savetxt(workdir + 'item_format_%d.dat' % (nstate),
            item_formats, fmt='%d')
    savetxt(workdir + 'item_names_%d.dat' % (nstate),
            item_names, fmt='%s')
    with open(workdir + 'item_def_doms_%d.dat' % (nstate), 'w') as f:
        for idd in item_def_doms:
            # if the parameter does not exist
            if idd == []:
                idd = [-2]
            for dd in idd:
                f.write(str(dd) + ' ')
            f.write('\n')

    # a = search_block(f, TAGS, 'ROOT', verbose=1, inv_TAGS=inv_TAGS)

    f.close()

    return ndfiles, elfiles
