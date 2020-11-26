/*  
 *  Copyright 2020 Francesco Gozzini < gozzini AT cpt.univ-mrs.fr >
 *
 *  This file is part of SL2CFOAM-NEXT.
 *
 *  SL2CFOAM-NEXT is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  SL2CFOAM-NEXT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with SL2CFOAM-NEXT. If not, see <http://www.gnu.org/licenses/>.
 */

#include "sl2cfoam_tensors.h"

#include "utils.h"
#include "common.h"

#define check(name, cond) \
    if (!(cond)) { fprintf(stderr, "test '%s' failed\n", name); exit(EXIT_FAILURE); }

// Testing a tensor with 4 indices, index n has size n^2 + 1
TENSOR_INIT(test, 4);

int main() {

    printf("Testing tensors...\n");

    tensor_ptr(test) t;
    TENSOR_CREATE(test, t, 4, 2, 5, 10, 17);

    double v;
    for (int i4 = 0; i4 < 17; i4++) {
    for (int i3 = 0; i3 < 10; i3++) {
    for (int i2 = 0; i2 < 5; i2++) {
    for (int i1 = 0; i1 < 2; i1++) {

        v = (double)((i1+1) * (i2+1) * (i3+1) * (i4+1));
        TENSOR_SET(v, t, 4, i1, i2, i3, i4);

    }
    }
    }
    }

    char tag[3][16];
    strcpy(tag[0], "sl2cfoam");
    strcpy(tag[1], "is");
    strcpy(tag[2], "awesome");
    TENSOR_TAG(t, tag, sizeof(tag));

    // store tensor
    TENSOR_STORE(t, "test/test_data/sl2cfoam_tensor_test.sl2t");

    tensor_ptr(test) tl;

    // load tensor
    TENSOR_LOAD(test, tl, 4, "test/test_data/sl2cfoam_tensor_test.sl2t");

    check("dim1", tl->dims[0] == 2);
    check("dim2", tl->dims[1] == 5);
    check("dim3", tl->dims[2] == 10);
    check("dim4", tl->dims[3] == 17);

    char decl[48];
    char (*tag_cast)[16] = (char (*)[16])t->tag;
    snprintf(decl, 48, "%s %s %s", tag_cast[0], tag_cast[1], tag_cast[2]);
    
    check("tag", strcmp(decl, "sl2cfoam is awesome") == 0);

    for (int i4 = 0; i4 < 17; i4++) {
    for (int i3 = 0; i3 < 10; i3++) {
    for (int i2 = 0; i2 < 5; i2++) {
    for (int i1 = 0; i1 < 2; i1++) {

        check("value get", TENSOR_GET(tl, 4, i1, i2, i3, i4) == TENSOR_GET(t, 4, i1, i2, i3, i4));

    }
    }
    }
    }

    TENSOR_FREE(t);
    TENSOR_FREE(tl);

    printf("... test passed.\n");

    // remove test file
    remove("test/test_data/sl2cfoam_tensor_test.sl2t");

    return EXIT_SUCCESS;
    
}