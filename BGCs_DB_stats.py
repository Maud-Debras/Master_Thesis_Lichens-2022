import csv
import xlsxwriter
import re

# Create the map from type to clusters
clusters_map_filename = "list-classes-rules.txt"
type_clusters_map = {}
clusters = []
with open(clusters_map_filename) as map_file:
    for line in map_file:
        cluster, types = line[:-1].split("\t")
        clusters.append(cluster)
        for t in types.split(" "):
            if not type_clusters_map.get(t):
                type_clusters_map[t] = []
            type_clusters_map[t].append(cluster)
cluster_index = {}
for cnt, cluster in enumerate(clusters):
    cluster_index[cluster] = cnt

# à changer avec un with: 
tsv_file = open("Clusters_all.tsv")
read_tsv = csv.reader(tsv_file, delimiter="\t")
regions_data = {}
for _, _, t, _, _, _, _, _, _, _, raw_region_name in read_tsv:
    region = raw_region_name[5:8]
    for cluster in type_clusters_map[t]:
        if not regions_data.get(region):
            regions_data[region] = [0] * len(clusters)
        regions_data[region][cluster_index[cluster]] += 1

workbook = xlsxwriter.Workbook('stats.xlsx')
worksheet = workbook.add_worksheet()

for cnt_cluster, cluster in enumerate(clusters):
    worksheet.write(0, 2+2*cnt_cluster, cluster)
    worksheet.write(0, 2+2*cnt_cluster+1, f"Abundance of {cluster}")

worksheet.write(0, 2+2*len(clusters)+1, "Total")

row = 2

with open("list-ids.txt") as f:
    for line in f:
        regions_group = line[:-1].split(",")
        for region in regions_group:
            worksheet.write(row, 0, region)
            region_values = regions_data[region]
            region_values_total = sum(region_values)
            worksheet.write(row, 2+2*len(clusters)+1, region_values_total)
            for cnt_val, val in enumerate(region_values):
                worksheet.write(row, 2+2*cnt_val, val)
                worksheet.write(row, 2+2*cnt_val+1, val/region_values_total)

            row += 1

workbook.close()
