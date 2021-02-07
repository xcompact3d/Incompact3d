/**********************************************************************************
 *
 * This file is part of Xcompact3d.
 *
 * Xcompact3d
 * Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
 * eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
 * 
 *     Xcompact3d is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation.
 * 
 *     Xcompact3d is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with the code.  If not, see <http://www.gnu.org/licenses/>.
 * -------------------------------------------------------------------------------
 * -------------------------------------------------------------------------------
 *     We kindly request that you cite Xcompact3d/Incompact3d in your
 *     publications and presentations. The following citations are suggested:
 * 
 *     1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
 *     incompressible flows: a simple and efficient method with the quasi-spectral
 *     accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
 * 
 *     2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
 *     problems with up to 0(10^5) computational cores, Int. J. of Numerical
 *     Methods in Fluids, vol 67 (11), pp 1735-1757
 *
 *********************************************************************************/

/* 
 * SLURM can send the signal SIGUSR2 180 seconds before the end of the simulation
 *
 * #SBATCH --signal=SIGUSR2@180
 *
 */

#include <signal.h>

typedef void (*sighandler_t)(int);

void signal_( int* signum, sighandler_t handler)
{
	signal(*signum, handler);
}
