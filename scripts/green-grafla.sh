#!/bin/sh

mkdir grap-files

export2graphlan.py --skip_rows 1,2 -i mpa-files/combine.mpa --tree grap-files/merged_abundance.tree.txt \
--annotation grap-files/merged_abundance.annot.txt --most_abundant 150 --annotations 2 \
--external_annotations 6 --abundance_threshold 15 \
--annotation_legend_font_size 8 --def_font_size 30
echo Output files saved inside grap-files folder

echo Color for Bacteroidetes changed from $#2d19ff to $#e6ab02
sed 's/#2d19ff/#e6ab02/g' grap-files/merged_abundance.annot.txt > temp.txt && mv temp.txt grap-files/merged_abundance.annot.txt

echo Color for Actinobacteria changed from $#29cc36 to $#e7298a
sed 's/#29cc36/#1b9e77/g' grap-files/merged_abundance.annot.txt > temp.txt && mv temp.txt grap-files/merged_abundance.annot.txt

echo Color for Acidobacteria changed from $#ff3333 to $#d95f03
sed 's/#ff3333/#e7298a/g' grap-files/merged_abundance.annot.txt > temp.txt && mv temp.txt grap-files/merged_abundance.annot.txt

echo Color for Firmicutes changed from $#00bfff to $#1b9e77
sed 's/#00bfff/#d95f03/g' grap-files/merged_abundance.annot.txt > temp.txt && mv temp.txt grap-files/merged_abundance.annot.txt

echo Color for Proteobacteria changed from $#00ff80 to $#7570b3
sed 's/#ffea00/#7570b3/g' grap-files/merged_abundance.annot.txt > temp.txt && mv temp.txt grap-files/merged_abundance.annot.txt


graphlan_annotate.py --annot grap-files/merged_abundance.annot.txt grap-files/merged_abundance.tree.txt grap-files/merged_abundance.xml

echo Generating the .png file
graphlan.py --dpi 300 --size 10 grap-files/merged_abundance.xml green-graphlan_graph.png --external_legends