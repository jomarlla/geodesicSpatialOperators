--Tested with:
--Oracle Database 12c Enterprise Edition Release 12.2.0.1.0 - 64bit Productionsystems
--Oracle Database 18c Enterprise Edition Release 18.0.0.0.0 - Production Version 18.3.0.0.0

--You must open a connection to any database and run the followings SQL commands:


/* --------------------------------------------------------------------------------------------------------
  Checking point to line minimum distance (the two first cases are AB2 and AB3 from the previous example)
  -------------------------------------------------------------------------------------------------------- 
  
  https://docs.oracle.com/database/121/SPATL/sdo_geom-sdo_distance.htm#SPATL1117
  */
CREATE TABLE LINES2 (NAME VARCHAR(3) PRIMARY KEY, GEOM SDO_GEOMETRY);
INSERT INTO LINES2 VALUES('AB1', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(5,52,6,51.4) ) );
INSERT INTO LINES2 VALUES('AB2', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(29,42,-77,39) ) );
INSERT INTO LINES2 VALUES('AB3', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(29,42,-70,-35) ) );


CREATE TABLE POINTS (
  NAME VARCHAR(3) PRIMARY KEY,
  GEOM SDO_GEOMETRY);
  
INSERT INTO USER_SDO_GEOM_METADATA (TABLE_NAME, COLUMN_NAME, DIMINFO, SRID) 
    VALUES ('POINTS','GEOM', 
       MDSYS.SDO_DIM_ARRAY 
        (MDSYS.SDO_DIM_ELEMENT('Longitude', -180, 180, 0.05), --Tolerance = 0.05 (*)
         MDSYS.SDO_DIM_ELEMENT('Latitude', -90, 90, 0.05) ),  --Tolerance = 0.05 (*)
       NULL);

INSERT INTO POINTS VALUES('X1', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,1,1), SDO_ORDINATE_ARRAY(5.5, 52) ) );
INSERT INTO POINTS VALUES('X2', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,1,1), SDO_ORDINATE_ARRAY(-22, 64) ) );
INSERT INTO POINTS VALUES('X3', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,1,1), SDO_ORDINATE_ARRAY(-22, 64) ) );

SELECT SDO_GEOM.SDO_DISTANCE(l.geom, p.geom,0.001) FROM lines2 l, points p WHERE l.name = 'AB1' and p.name = 'X1';
SELECT SDO_GEOM.SDO_DISTANCE(l.geom, p.geom,0.001) FROM lines2 l, points p WHERE l.name = 'AB2' and p.name = 'X2';
SELECT SDO_GEOM.SDO_DISTANCE(l.geom, p.geom,0.001) FROM lines2 l, points p WHERE l.name = 'AB3' and p.name = 'X3';

--Case 1: 23723,6561485326
--Case 2: 1013727,29054293
--Case 3: 3923508,12609771
