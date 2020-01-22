/*
 * Experiments for the paper [1]
 * [1] New paper (sending status)
 *
 * Copyright (c) Jose Martinez-Llario (2018-2019) and licensed
 * under the GPL-3.0-or-later License. 
 * The algorithm uses the library geographiclib from Charles Karney to
 * perform the direct and inverse geodetic problems
 * 
 *  This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,e
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 **********************************************************************/
	
/*
   Before running this file run in PostGIS the file "table_geodesicintersection_several_software.sql"
   and read the instructions there.
*/

-- Print out all the X coordinates in DMS and decimal degrees 
 
		/* -- S.12345 (0.3mm in equador) S.12345678(0,3um in equador)
		SELECT i.ncase, i.software, st_astext(geom) as degrees, (ST_AsLatLonText(geom, 'D°M''S.SSSSSSSS"')) as DMS from inter i order by i.software, i.ncase;	
				

		"ncase","software","degrees","dms"
		"Case 1","ArcGIS","POINT(14.5286278725001 54.0572557445001)","54°3'26.12068020"" 14°31'43.06034100"""
		"Case 2","ArcGIS","POINT(5.22744828600003 51.865640108)","51°51'56.30438880"" 5°13'38.81382960"""
		"Case 3","ArcGIS","POINT(-14.563892229 54.7170126655001)","54°43'1.24559580"" -14°33'50.01202440"""
		"Case 4","ArcGIS","POINT(-79.2828199255 50.4790691370001)","50°28'44.64889320"" -79°16'58.15173180"""
		"Case 1","ArcGISUTM","POINT(14.5285483925224 54.0573011622007)","54°3'26.28418392"" 14°31'42.77421308"""
		"Case 2","ArcGISUTM","POINT(5.22744847138889 51.8656407136111)","51°51'56.30656900"" 5°13'38.81449700"""
		"Case 3","ArcGISUTM","POINT(-14.5638559686111 54.7170296177778)","54°43'1.30662400"" -14°33'49.88148700"""
		"Case 4","ArcGISUTM","POINT(-79.2828200318611 50.4790696119444)","50°28'44.65060300"" -79°16'58.15211470"""
		"Case 1","Oracle&Google","POINT(14.5285478358549 54.0573010333795)","54°3'26.28372017"" 14°31'42.77220908"""
		"Case 2","Oracle&Google","POINT(5.22745698883988 51.8656617779676)","51°51'56.38240068"" 5°13'38.84515982"""
		"Case 3","Oracle&Google","POINT(-14.544131677975 54.6717159071832)","54°40'18.17726586"" -14°32'38.87404071"""
		"Case 4","Oracle&Google","POINT(-79.2501486314468 50.4442144186783)","50°26'39.17190724"" -79°15'0.53507321"""
		"Case 1","PG_Paper","POINT(14.5285478498717 54.0573013091912)","54°3'26.28471309"" 14°31'42.77225954"""
		"Case 2","PG_Paper","POINT(5.22745711452157 51.8656654013764)","51°51'56.39544495"" 5°13'38.84561228"""
		"Case 3","PG_Paper","POINT(-14.5638557443078 54.7170296089477)","54°43'1.30659221"" -14°33'49.88067951"""
		"Case 4","PG_Paper","POINT(-79.2828016866239 50.4790974467668)","50°28'44.75080836"" -79°16'58.08607185"""


		*/
	
								  
									  
--TestB
 create or replace view distancias as
 select l.name as line, i.ncase, i.software,
 st_distance (
    st_project (st_startpoint (l.geom)::geography,
      st_distance(st_startpoint (l.geom)::geography, i.geom::geography),
	  st_azimuth(st_startpoint (l.geom)::geography, st_endpoint (l.geom)::geography))
    ,
	i.geom) as distancia
   from lines l, inter  i where l.ncase = i.ncase order by i.software, l.name, l.ncase;
   
   /*
     The distances for the Paper algorithm are ecaxtly 0, what means we can not check with double precision. We will use
	 the software Maxima with exact computation and 60 digits of precision to check the results.
   */
	/* select * from distancias;

			"line","ncase","software","distancia"
			"AB","Case 1","ArcGIS","6.45663129"
			"AB","Case 2","ArcGIS","2.4525944"
			"AB","Case 3","ArcGIS","2.11612342"
			"AB","Case 4","ArcGIS","0.78881627"
			"CD","Case 1","ArcGIS","1.83208025"
			"CD","Case 2","ArcGIS","1.81710257"
			"CD","Case 3","ArcGIS","2.83863448"
			"CD","Case 4","ArcGIS","3.09717648"
			"AB","Case 1","ArcGISUTM","0.03870534"
			"AB","Case 2","ArcGISUTM","2.39510984"
			"AB","Case 3","ArcGISUTM","0.00049673"
			"AB","Case 4","ArcGISUTM","0.75182667"
			"CD","Case 1","ArcGISUTM","0.00326752"
			"CD","Case 2","ArcGISUTM","1.77247583"
			"CD","Case 3","ArcGISUTM","0.01333927"
			"CD","Case 4","ArcGISUTM","3.04409387"
			"AB","Case 1","Oracle&Google","0.00776586"
			"AB","Case 2","Oracle&Google","0.29705727"
			"AB","Case 3","Oracle&Google","4887.9260394"
			"AB","Case 4","Oracle&Google","4148.96691918"
			"CD","Case 1","Oracle&Google","0.02716831"
			"CD","Case 2","Oracle&Google","0.30933963"
			"CD","Case 3","Oracle&Google","446.46370698"
			"CD","Case 4","Oracle&Google","3965.88536726"
			"AB","Case 1","PG_Paper","0"
			"AB","Case 2","PG_Paper","0"
			"AB","Case 3","PG_Paper","0"
			"AB","Case 4","PG_Paper","0"
			"CD","Case 1","PG_Paper","0"
			"CD","Case 2","PG_Paper","0"
			"CD","Case 3","PG_Paper","0"
			"CD","Case 4","PG_Paper","0"
			

	*/   
	
	select d.software, d.ncase, sqrt ((var_pop (d.distancia) + avg(d.distancia)*avg(d.distancia))*count(*)) as distanciacuadratica,
	       max(d.distancia) as distanciamax
	   from distancias d group by d.software, d.ncase order by software, ncase;
			 /*
				"software","ncase","distanciacuadratica","distanciamax"
				"ArcGIS","Case 1","6.71152781842161","6.45663129"
				"ArcGIS","Case 2","3.05238939862069","2.4525944"
				"ArcGIS","Case 3","3.54059656551228","2.83863448"
				"ArcGIS","Case 4","3.19604963291904","3.09717648"
				"ArcGISUTM","Case 1","0.0388430177955575","0.03870534"
				"ArcGISUTM","Case 2","2.97963449328924","2.39510984"
				"ArcGISUTM","Case 3","0.0133485154540046","0.01333927"
				"ArcGISUTM","Case 4","3.13556228308335","3.04409387"
				"Oracle&Google","Case 1","0.0282564266989954","0.02716831"
				"Oracle&Google","Case 2","0.428875306293554","0.30933963"
				"Oracle&Google","Case 3","4908.27370959433","4887.9260394"
				"Oracle&Google","Case 4","5739.52726648258","4148.96691918"
				"PG_Paper","Case 1","0","0"
				"PG_Paper","Case 2","0","0"
				"PG_Paper","Case 3","0","0"
				"PG_Paper","Case 4","0","0"


			*/

 --Test A
 create view angulos as
  select l.name as line, i.ncase, i.software,
    ( st_azimuth(st_startpoint (l.geom)::geography, i.geom::geography) 
      - st_azimuth (st_startpoint (l.geom)::geography, st_endpoint (l.geom)::geography)) * (180 / pi()) * 60 * 60
	  as azimuth

 from lines l, inter i where l.ncase = i.ncase order by i.software, l.name, l.ncase;

		/*
		 select * from angulos;
		 
		"line","ncase","software","azimuth"
		"AB","Case 1","ArcGIS","200.470010592739"
		"AB","Case 2","ArcGIS","23.3788048870166"
		"AB","Case 3","ArcGIS","-0.13273993637841"
		"AB","Case 4","ArcGIS","0.0825470657176539"
		"CD","Case 1","ArcGIS","39.0777507459199"
		"CD","Case 2","ArcGIS","5.792965368055"
		"CD","Case 3","ArcGIS","-0.120008371313372"
		"CD","Case 4","ArcGIS","-0.102597571792648"
		"AB","Case 1","ArcGISUTM","1.20113893692031"
		"AB","Case 2","ArcGISUTM","22.830885709543"
		"AB","Case 3","ArcGISUTM","-3.11587161368441e-05"
		"AB","Case 4","ArcGISUTM","0.0786762229322035"
		"CD","Case 1","ArcGISUTM","-0.069745949324155"
		"CD","Case 2","ArcGISUTM","5.65068972716473"
		"CD","Case 3","ArcGISUTM","-0.000563941479540475"
		"CD","Case 4","ArcGISUTM","-0.100839148441605"
		"AB","Case 1","Oracle&Google","0.240998039588809"
		"AB","Case 2","Oracle&Google","2.83179332312147"
		"AB","Case 3","Oracle&Google","-306.751586967694"
		"AB","Case 4","Oracle&Google","434.550652428659"
		"CD","Case 1","Oracle&Google","0.579914835069217"
		"CD","Case 2","Oracle&Google","0.986151863827824"
		"CD","Case 3","Oracle&Google","-18.8879668779382"
		"CD","Case 4","Oracle&Google","-131.36525892841"
		"AB","Case 1","PG_Paper","-4.24794883256267e-09"
		"AB","Case 2","PG_Paper","-7.69439788539654e-09"
		"AB","Case 3","PG_Paper","2.28999937065373e-11"
		"AB","Case 4","PG_Paper","8.01499779728806e-11"
		"CD","Case 1","PG_Paper","-1.39231961735747e-08"
		"CD","Case 2","PG_Paper","-3.16019913150215e-09"
		"CD","Case 3","PG_Paper","-1.7174995279903e-11"
		"CD","Case 4","PG_Paper","6.8699981119612e-11"
		*/

 select d.software, d.ncase, max(d.azimuth) as azimuthmax
	   from angulos d group by d.software, d.ncase order by software, ncase;

	 /*
		"software","ncase","azimuthmax"
		"ArcGIS","Case 1","200.470010592739"
		"ArcGIS","Case 2","23.3788048870166"
		"ArcGIS","Case 3","-0.120008371313372"
		"ArcGIS","Case 4","0.0825470657176539"
		"ArcGISUTM","Case 1","1.20113893692031"
		"ArcGISUTM","Case 2","22.830885709543"
		"ArcGISUTM","Case 3","-3.11587161368441e-05"
		"ArcGISUTM","Case 4","0.0786762229322035"
		"Oracle&Google","Case 1","0.579914835069217"
		"Oracle&Google","Case 2","2.83179332312147"
		"Oracle&Google","Case 3","-18.8879668779382"
		"Oracle&Google","Case 4","434.550652428659"
		"PG_Paper","Case 1","-4.24794883256267e-09"
		"PG_Paper","Case 2","-3.16019913150215e-09"
		"PG_Paper","Case 3","2.28999937065373e-11"
		"PG_Paper","Case 4","8.01499779728806e-11"

	*/

 


-- Geodesic distance to the exact coordinates from the proposed algorithm
	select i1.ncase, i1.software, st_distance (i0.geom::geography, i1.geom::geography) 
		from inter i0, inter i1 where i0.software = 'PG_Paper' and i1.software <> 'PG_Paper' and
								  i0.ncase = i1.ncase order by i1.software, i1.ncase;
								  
				/*
					"ncase","software","st_distance"
					"Case 1","ArcGIS","7.29263703"
					"Case 2","ArcGIS","2.87921965"
					"Case 3","ArcGIS","3.01423603"
					"Case 4","ArcGIS","3.40485429"
					"Case 1","ArcGISUTM","0.03912131"
					"Case 2","ArcGISUTM","2.81065993"
					"Case 3","ArcGISUTM","0.01448832"
					"Case 4","ArcGISUTM","3.35896991"
					"Case 1","Oracle&Google","0.03071324"
					"Case 2","Oracle&Google","0.40325079"
					"Case 3","Oracle&Google","5202.08509206"
					"Case 4","Oracle&Google","4520.25667908"
				*/			
				
--Geodetic Distance beetween start points and intersection point
	select i.software, l.name, st_distance(st_startpoint(l.geom), i.geom::geography) as AX_CX 
		from lines l, inter i 
		   where i.software = 'PG_Paper' and l.ncase = i.ncase order by i.software, l.name, i.ncase;
	/*
		"software","name","ax_cx"
		"PG_Paper","AB","6646.65565016"
		"PG_Paper","AB","21637.10319232"
		"PG_Paper","AB","3454490.28172529"
		"PG_Paper","AB","2003881.86798765"
		"PG_Paper","CD","9663.2474396"
		"PG_Paper","CD","64703.2463332"
		"PG_Paper","CD","5558129.37316681"
		"PG_Paper","CD","11347603.1157086"
	*/
								  
								  
								  