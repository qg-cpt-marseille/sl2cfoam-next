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

#ifndef __SL2CFOAM_FOLDERS_H__
#define __SL2CFOAM_FOLDERS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

////////////////////////////////////////////////////////////////////////
//  Manage folders configuration.
////////////////////////////////////////////////////////////////////////

// Contains the folders structure.
struct sl2cfoam_data_folders {

    char* vertex;
    char* vertex_imm;
    char* vertex_imm_boost;
    char* vertex_imm_ampl;

};

// Returns the folder structure.
struct sl2cfoam_data_folders* sl2cfoam_get_data_folders();

// Call to free memory of subfolders paths.
void sl2cfoam_free_data_folders(struct sl2cfoam_data_folders* df);

// Check if folders exist/can be created.
void sl2cfoam_check_data_folders(struct sl2cfoam_data_folders* df);


#ifdef SL2CFOAM_SHORT_NAMES
#define get_data_folders(...) sl2cfoam_get_data_folders(__VA_ARGS__)
#define free_data_folders(...) sl2cfoam_free_data_folders(__VA_ARGS__)
#define check_data_folders(...) sl2cfoam_check_data_folders(__VA_ARGS__)
#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_FOLDERS_H__*/
