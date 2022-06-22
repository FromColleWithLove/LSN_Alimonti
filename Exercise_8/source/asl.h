/*
	asl : Alimonti Statistics Library
	This library contains methods to perform needed calculation for LSN course.
	Copyright (C) 2022  Alimonti Davide
	
	This program is free software: you can redistribute it and/or modify
    	it under the terms of the GNU General Public License as published by
    	the Free Software Foundation, either version 3 of the License, or
    	(at your option) any later version.

    	This program is distributed in the hope that it will be useful,
    	but WITHOUT ANY WARRANTY; without even the implied warranty of
    	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    	GNU General Public License for more details.

    	You should have received a copy of the GNU General Public License
    	along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

extern void test();
double arraysum(int d,double a[]);
double arraysum2(int d,double a[]);
extern double avg(int d,double a[]);
extern double avg2(int d,double a[]);
extern double sq_mod(int d,double a[]);
extern double pearson_rho(int d,double a[],double b[]);
extern double autocorr_disc(double a[],int tmax, int t);
extern double Error(double sum, double sum2, int iblk);
