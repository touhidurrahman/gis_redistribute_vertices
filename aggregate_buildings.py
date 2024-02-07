##USER INPUTS
    
alkis_path = r"E:\OneDrive - VEI\Spekter\aggregate_buildings\aggregateBuildings_sample\buildings.shp"
working_directory = r"E:\OneDrive - VEI\Spekter\aggregate_buildings"
user_distance = 1  #CHANGE HERE IF NEEDED
aggregation_distance = "0.3 Meters"


##PLEASE DO NOT TOUCH THE FOLLOWING CODE IF NOT NEEDED

##PYTHON LIBRARIES

import geopandas as gpd
from shapely.geometry import Point, MultiPoint
from shapely.geometry import LineString
import arcpy
import os
import shutil
import time
import sys


user_distance = float(user_distance)


start_time = time.time()

##STORING TEMPORARY FILES
arcpy.env.overwriteOutput = True

try:
    print("Creating temporary folders...")

    output = f"{working_directory}\\output_folder"
    if not os.path.exists(output):
        os.makedirs(output)

    temp1 = f"{working_directory}\\temp1"
    if not os.path.exists(temp1):
        os.makedirs(temp1)

    temp2 = f"{working_directory}\\temp2"
    if not os.path.exists(temp2):
        os.makedirs(temp2)

    temp3 = f"{working_directory}\\temp3"
    if not os.path.exists(temp3):
        os.makedirs(temp3)
except Exception as e:
    print(f"Failed to create temporary folders: {str(e)}")


##AGGREGATE BUILDINGS
try:
    print("Aggregating buildings...")
    arcpy.cartography.AggregatePolygons(
        in_features= f"{alkis_path}",
        out_feature_class= f"{output}\\buildings_aggregate.shp",
        aggregation_distance= f"{aggregation_distance}",
        minimum_area="0 SquareMeters",
        minimum_hole_size="0 SquareMeters",
        orthogonality_option="NON_ORTHOGONAL",
        barrier_features=None,
        aggregate_field=None
    )
    buildings_aggregate_path = f"{output}\\buildings_aggregate.shp"
except Exception as e:
    print(f"Failed to aggregate buildings: {str(e)}")
    

##SIMPLIFY POLYGONS
    
try:
    print("Simplifying aggregated buildings...")
    arcpy.cartography.SimplifyPolygon(
        in_features=f"{buildings_aggregate_path}",
        out_feature_class=f"{output}\\buildings_simplify.shp",
        algorithm="POINT_REMOVE",
        tolerance="0.2 Meters",
        minimum_area="0 SquareMeters",
        error_option="RESOLVE_ERRORS",
        collapsed_point_option="NO_KEEP",
        in_barriers=None
    )

    buildings_simplify_path = f"{output}\\buildings_simplify.shp"
except Exception as e:
    print(f"Failed to simplify aggregated buildings: {str(e)}")


##REDISTRIBUTING VERTICES
    
try:
    print("Initializing vertice redistribution algorithm...")
    gdf = gpd.read_file(buildings_simplify_path)

    k = 1
    for index, row in gdf.iterrows():
        print(f"Progress: {k}")
        polygon = row['geometry']
        vertices = list(polygon.exterior.coords)
        
        all_coordinates = ""
        for i in range (len(vertices)-1):
            vertice1 = vertices[i]
            vertice2 = vertices[i+1]
            line = LineString([vertice1, vertice2])
            line_length = line.length
            print(f"Line Length: {line_length}")
            node = round(line_length/user_distance)
            if node == 0:
                node = 1
            if node>0:
                print(f"Node Number: {node}")
                final_distance = line_length/node

                new_coordinates = [line.interpolate(i * final_distance).coords[0] for i in range(1, node)]

                if i == 0:
                    all_coordinates = [vertice1] + new_coordinates
                else:
                    if all_coordinates == "":
                        all_coordinates = [vertice1] + new_coordinates
                    else:
                        all_coordinates = all_coordinates + [vertice1] + new_coordinates
        
        # Create Point objects for each coordinate
        points = [Point(coords) for coords in all_coordinates]
        
        print(points)

        # Create a GeoDataFrame with the Point geometries
        gdf = gpd.GeoDataFrame(geometry=points)
        
        # Save the GeoDataFrame to a shapefile
        output_shapefile = f"{temp1}\\point_{k}.shp"
        gdf.to_file(output_shapefile)
        
        k = k + 1
        print("---------------")
except Exception as e:
    print(f"Failed to redistribute vertices: {str(e)}")
    print(f"Failed to create Point_{k}.shp")
    sys.exit(1)

##POINT TO LINE
try:
    print("Creating lines from points...")
    arcpy.env.workspace = temp1
    pointFeatures = arcpy.ListFeatureClasses()
    pointFeatures_count = len(pointFeatures)

    j = 1
    for point in pointFeatures:
        arcpy.management.PointsToLine(
            Input_Features=f"{point}",
            Output_Feature_Class=f"{temp2}//line_{j}.shp",
            Line_Field=None,
            Sort_Field=None,
            Close_Line="CLOSE",
            Line_Construction_Method="CONTINUOUS",
            Attribute_Source="NONE",
            Transfer_Fields=None
        )
        print(f"Progress: {j} of {pointFeatures_count}")
        j = j + 1
        print(f"Line generated for {point}")
except Exception as e:
    print(f"Failed to create line from point: {str(e)}")
    print(f"Failed to create line from point at progress: {j}")
    sys.exit(1)

##LINE TO POLYGON
try:
    print("Creating polygons from lines...")
    arcpy.env.workspace = temp2
    lineFeatures = arcpy.ListFeatureClasses()
    lineFeatures_count = len(lineFeatures)

    l = 1
    for line in lineFeatures:
        arcpy.management.FeatureToPolygon(
            in_features= f"{line}",
            out_feature_class=f"{temp3}\\building_{l}.shp",
            cluster_tolerance=None,
            attributes="ATTRIBUTES",
            label_features=None
        )
        print(f"Progress: {l} of {lineFeatures_count}")
        l = l + 1
        print(f"Created building from {line}")
except Exception as e:
    print(f"Failed to create polygon from line: {str(e)}")
    print(f"Failed to create polygon from lines at progress: {l}")
    sys.exit(1)

##MERGE BUILDINGS
try:
    print("Merging buildings...")
    arcpy.env.workspace = temp3
    polygonFeatures = arcpy.ListFeatureClasses()

    polygonFeatures_input = '; '.join(polygonFeatures)

    arcpy.management.Merge(
        inputs=f"{polygonFeatures_input}",
        output=f"{output}\\final_buildings.shp",
        add_source="NO_SOURCE_INFO"
    )

    buildings_vertices_redistributed = f"{output}\\final_buildings.shp"
except:
    print(f"Failed to merge buildings: {str(e)}")
    sys.exit(1)


##REFERENCE THE FINAL FILE

try:
    print("Geo-referencing merged buildings...")
    spatial_ref = arcpy.Describe(alkis_path).spatialReference
    arcpy.management.DefineProjection(buildings_vertices_redistributed, spatial_ref)
except Exception as e:
    print(f"Failed to Geo-reference the merged buildings: {str(e)}")
    sys.exit(1)

##DELETE TEMPORARY FOLDERS
    
try:
    print("Deleting temporary folders...")
    if os.path.exists(temp1):
        shutil.rmtree(temp1)
    if os.path.exists(temp2):
        shutil.rmtree(temp2)
    if os.path.exists(temp3):
        shutil.rmtree(temp3)
except Exception as e:
    print(f"Failed to delete temporary folders: {str(e)}")

print("------------------")
print("PROCESS COMPLETED")
print(f"Check this folder: {output}")

##TIME ROCORDING
print(f"Time elapsed: {(time.time()-start_time)/60} minutes")