#don't use commas in snp names
import argparse
import time
import shutil
import glob
import os, sys
import pandas as pd
import math
import numpy as np
import subprocess
import warnings
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.interpolate import Rbf
import geopandas as gpd
from shapely.geometry import Point
from matplotlib.colors import LinearSegmentedColormap
import geopandas


#User input paths
parser = argparse.ArgumentParser()

parser.add_argument("-input_format", choices=["vcf","csv"],help="Option 'vcf' is used when the available input format are .vcf or vcf.gz files. The vcf files need to be quality filtered already! The 'csv' option is used for tab-delimited csv file of a dataframe with individuals as columns and variants as rows. "
  "Insert variant data as 'A' for ancestral, 'D' for derived allele, and 'X' for missing data. Add variant names as first column and sequence name as header row. "
  "For both options: Refrain from using comma separated marker or sequence names.")

parser.add_argument("-input", required=True, help="For 'vcf' give the folder to the vcf files. For 'csv', give the path to one csv file. If -freqmap is specified, this needs to specify phyloimputed (or external) tab-separated csv file for allele frequency map.")
parser.add_argument("-output", required=True, help="Path to output folder.")

parser.add_argument('-freqmap', action='store_true', help='Generate allele frequency map for specified SNP.')
parser.add_argument("-f_snp", help="Define SNP name to generate allele frequency maps.")
parser.add_argument("-f_coordinates", help="Provide path to csv file with sample names in first column (case-sensitive to input file) and coordinates in second column.")
parser.add_argument("-color", choices=["blue","orange","pink","red","green","yellow","purple","violet","grey"],help="Select color [default:blue]") 
parser.add_argument("-contour", type=str, help="Define resolution of allele frequencies. Default: 15")
parser.add_argument("-smoothing", type=str, help="Define interpolated missing data points between sampling coordinates. The higher the value, the more extreme the smoothing. Default: 2.3")
parser.add_argument("-ancestral_coordinates", action="store_true",help="Specify if marking coordinates with ancestral alleles of SNP in black circles.")
parser.add_argument("-derived_coordinates", action="store_true",help="Specify if marking coordinates with ancestral alleles of SNP in white circles.")

parser.add_argument("-whole_world", action="store_true",help="Specify if visualizing map of whole world.")
parser.add_argument("-continent",choices=["Oceania","Africa","North America","Asia","South America","Europe"],nargs='+',help="Select one or more continents that you want to present")
parser.add_argument("-country",nargs='+',help="Select one or more countries that you want to present. Nomenclature must align with the countries' nomenclature in file 'Countries_list.csv'.")
parser.add_argument("-af_map", default="svg", choices=["pdf","png","svg"],help="Specify output format of allele frequency map. Default: svg")


args = parser.parse_args()

current_dir = os.path.dirname(os.path.abspath(__file__))

print()
print(f"Start generating allele frequency map for: {args.f_snp}")
#open color maps
colors_df= pd.read_csv("color_maps.csv", sep="\t")
color_str = colors_df.loc[colors_df['color'] == args.color, 'color_map'].values[0]
# Convert string to list
color_list = [c.strip() for c in color_str.split(',')]
# Recreate the colormap
custom_cmap = LinearSegmentedColormap.from_list("custom", color_list)


coordinates= pd.read_csv(args.f_coordinates, sep="\t")
coordinates.columns=["sample","lat","lon"]


df= pd.read_csv(args.input, sep="\t")
df.set_index(df.columns[0], inplace=True)
df_snp = df.loc[[args.f_snp]]
df_snp = df_snp.T

merged_df = pd.merge(df_snp, coordinates, left_index=True, right_on='sample', how="inner")
merged_df.reset_index(drop=True, inplace=True)
merged_df["lat"] = merged_df["lat"].astype(float).round(1).astype(int)
merged_df["lon"] = merged_df["lon"].astype(float).round(1).astype(int)
merged_df[args.f_snp] = merged_df[args.f_snp].str.upper()

#exclude missing observations
merged_df = merged_df[merged_df[args.f_snp] != 'X']


merged_df = (
    merged_df.groupby([args.f_snp, 'lat', 'lon'])
      .size()
      .div(len(merged_df))
      .reset_index(name='weight')
)


#Prepare frequencies for derived alleles
merged_df_D = merged_df[merged_df[args.f_snp] == 'D']
merged_df_A = merged_df[merged_df[args.f_snp] == 'A']
# print(merged_df_D)

#2.1 PLOT
#Extracting x,y and values (z)
z=merged_df_D.weight
y=merged_df_D.lat
x=merged_df_D.lon

#Extract same info from ancestral alleles
z_A=merged_df_A.weight
y_A=merged_df_A.lat
x_A=merged_df_A.lon


# Custom weight function that decreases with distance
def custom_weight(distance):
    return np.exp(-distance**2)

# Extend the range of your regular grid
xi, yi = np.meshgrid(np.linspace(x.min() - 4, x.max() + 4, 100), np.linspace(y.min() - 4, y.max() + 4, 100))

z = gaussian_filter(z, sigma=5) #more or less smoothing
# Create a radial basis function (RBF) interpolator with custom weighting
rbf = Rbf(x, y, z, function="gaussian",epsilon=float(args.smoothing)) #2.3)


# Perform the RBF interpolation on the grid
zi = rbf(xi, yi)
# Ensure interpolated values are not higher than the original values
# zi = np.minimum(zi, np.max(z))
# zi = np.maximum(zi, np.min(z))
# Set negative values in zi to 0
zi[zi < 0] = 0


# Set the precision for the NumPy array to 3 decimal places
np.set_printoptions(precision=3, suppress=True)


# Create a GeoDataFrame with the interpolated values on the grid
grid_geometry = [Point(lon, lat) for lon, lat in zip(xi.flatten(), yi.flatten())]
gdf_interpolated = gpd.GeoDataFrame(geometry=grid_geometry, data={'Interpolated_Value': zi.flatten()})


# Load a world map from GeoPandas
if args.whole_world:
    world = gpd.read_file('./ne_110m_land/ne_110m_land.shp')
    s_A = 3
else:
    world = gpd.read_file('./ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')
    if args.continent:
        world = world[world["CONTINENT"].isin(args.continent)]
        # polygon = world[world["CONTINENT"].isin(args.continent)].union_all() #for making clean corners to ocean, flawed coordinate system
        s_A = 8
    elif args.country:
        world = world[world["SOVEREIGNT"].isin(args.country)]
        # polygon = world[world["SOVEREIGNT"].isin(args.country)].union_all() #for making clean corners to other countries, flawed coordinate system
        s_A = 10
    else:
        s_A = 3

#to keep marking within country border:
# if args.border == "free":
# if args.country:
    # mask = np.array([
        # [polygon.contains(Point(lon, lat)) for lon, lat in zip(row_xi, row_yi)]
        # for row_xi, row_yi in zip(xi, yi)
    # ])

    # Apply the mask: set values outside the polygon to np.nan
    # zi = np.where(mask, zi, np.nan)


# Plot the world map with the zoomed region
fig, ax = plt.subplots(figsize=(12, 8))


# print(world)



ax.set_aspect('equal')
# Plot borders
world.boundary.plot(ax=ax, color='black', linewidth=0.5)


# Plot the interpolated values on the world map with the custom colormap and transparency    

plt.xticks([])  # Remove x-axis ticks and labels
plt.yticks([])  # Remove y-axis ticks and labels
plt.title(args.f_snp)

plt.contourf(xi, yi, zi, levels=int(args.contour), cmap=custom_cmap, alpha=1, extend='both')
plt.colorbar(label='Relative SNP frequencies')

#Add ancestral alleles
if args.ancestral_coordinates == True:
    plt.scatter(x_A, y_A, c=z_A, marker='.', edgecolor='black', s=2*s_A, linewidth=(2*s_A)/10, label='ancestral alleles')
    legend = plt.legend(loc='upper right', frameon=True, title="Coordinates with")
    legend.get_frame().set_facecolor('#999999')

if args.derived_coordinates:
    plt.scatter(x, y, c=z, marker='.', s=s_A, edgecolor='white', linewidth=s_A/10, label='derived alleles')
    legend = plt.legend(loc='upper right', frameon=True, title="Coordinates with")
    legend.get_frame().set_facecolor('#999999')
    

if args.af_map == "svg":
    plt.savefig(str(args.output)+'_allelefreq_map_'+str(args.f_snp)+'.svg', format='svg',transparent=True)
elif args.af_map == "png":
    plt.savefig(str(args.output)+'_allelefreq_map_' + str(args.f_snp) + '.png', format='png', dpi=300, transparent=True)
elif args.af_map == "pdf":
    plt.savefig(str(args.output)+'_allelefreq_map_'+str(args.f_snp)+'.pdf', format='pdf',transparent=True)

print()
print("Finished generating allele frequency map!")
