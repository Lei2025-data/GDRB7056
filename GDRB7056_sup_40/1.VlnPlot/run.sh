/home/qrjia/anaconda3/bin/python /home/qrjia/exercise/PythonScript/conv.py --conv xlsx2txt --fr ../1.Violin/23120545增补基因.xlsx --to ../1.Violin/marker.glist
sed 1d ../1.Violin/marker.glist > ../1.Violin/marker.txt
/Bio/bin/Rscript-3.5.1_conda work.r
