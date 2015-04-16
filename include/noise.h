/*
 * Copyright (c) 2006 Filipe Maia
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


#ifndef _NOISE_H_
#define _NOISE_H_ 1
#ifdef GSL_FOUND
/* returns a random number from the poisson distribution of mean L */
void init_random_generator();
int get_poisson_random_number(double L);
int get_poisson_gaussian_approximate_random_number(double L);
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void generate_poisson_noise(CCD * det);
  void generate_gaussian_noise(CCD * det);
#ifdef __cplusplus
}
#endif

#endif
