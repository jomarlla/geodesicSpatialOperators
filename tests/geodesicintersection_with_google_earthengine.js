// Google earth enginee
// Last version tested: march 2019

// Run this code into the google online editor:
// https://code.earthengine.google.com/#

//WGS84 by default
var AB1 = ee.Geometry.LineString ( [[14.5, 54], [14.6, 54.2]] );
var CD1 = ee.Geometry.LineString ( [[14.4, 54.1], [14.7, 54]] );

var AB2 = ee.Geometry.LineString ( [[5, 52], [6, 51.4]] );
var CD2 = ee.Geometry.LineString ( [[4.5, 51.5], [5.5, 52]] );

var AB3 = ee.Geometry.LineString ( [[29, 42], [-77,39]] );
var CD3 = ee.Geometry.LineString ( [[0, 6], [-22,64]] );

var AB4 = ee.Geometry.LineString ( [[-92, 35], [52, 40]] );
var CD4 = ee.Geometry.LineString ( [[20, -8], [-95, 49]] );


var px1, px2, px3, px4;
px1 = AB1.intersection(CD1, ee.ErrorMargin(0.001));
px2 = AB2.intersection(CD2, ee.ErrorMargin(0.001));
px3 = AB3.intersection(CD3, ee.ErrorMargin(0.001));
px4 = AB3.intersection(CD3, ee.ErrorMargin(0.001));

print (px1.coordinates());
print (px2.coordinates());
print (px3.coordinates());
print (px4.coordinates());

//Case 1: 14.528547835854896,54.05730103337946
//Case 2: 5.227456988839927,51.865661777967624
//Case 3: 14.544131677975026,54.67171590718321
//Case 4: 79.25014863144679,50.44421441867829



