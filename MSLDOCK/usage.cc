/* MSLDOCK version 0.4 
 edited by Chao Li in 20200807
 Added some version information and general usage for MSLDOCK
 Following is original notes for Autodock
-----------------------------------------------------------------------
*/

/*

 $Id: usage.cc,v 1.6 2010/08/27 00:05:09 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful \
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "usage.h"


/*----------------------------------------------------------------------------*/
void usage( FILE *const file, const char *const programname )
/*----------------------------------------------------------------------------*/
{
    const char AutoDockHelp[] = \
    "\t-p parameter_filename\n" \
    "\t\t\t-l log_filename\n" \
    "\t\t\t-k (Keep original residue numbers)\n" \
    "\t\t\t-i (Ignore header-checking)\n" \
    "\t\t\t-t (Parse the PDBQT file to check torsions, then stop.)\n" \
    "\t\t\t-d (Increment debug level)\n" \
    "\t\t\t-C (Print copyright notice)\n" \
    "\t\t\t--version (Print msldock version)\n" \
    "\t\t\t--help (Display this message)\n" \
    "\n" \
    "\tThe new method added in MSLDOCK now only support command-line usage!\n" \
    "\tCan be used to run single docking tasks in parallel\n" \
    "\t\t(not used to parallelize multiple tasks!)\n" \
    "\tThe usage for other methods originally contained in AutoDock 4.2.6 is unchanged\n" \
    "\tSet run_mode in .dpf file to \"msldock_run\" to start docking with MSLRDPSO\n" \
    "\tGeneral usage for MSLDOCK:\n" \
    "\t\tNormal mode (Serial mode): msldock -p file.dpf -l file.dlg\n" \
    "\t\tParallel mode: msldock.omp -p file.dpf -l file.dlg\n\n"; 

    fprintf(file, "usage: %s %s\n", programname, AutoDockHelp);
}
/* EOF */
