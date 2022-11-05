import os,sys
# use like that:
# python dataset_paths.py name_of_txt_file T2_DE_RWTH

site_redirector = {
    'T2_DE_RWTH' : 'root://grid-cms-xrootd.physik.rwth-aachen.de:1094/',
    'T2_DE_DESY' : 'root://dcache-cms-xrootd.desy.de:1094/'
}

txt_name = "doublemuon2017"
site = "T2_DE_RWTW"

if len(sys.argv) > 1:
    txt_name = sys.argv[1]
if len(sys.argv) > 2:
    site = sys.argv[2]
    
fIN = open("%s.txt" % (txt_name),'r')
fOUT = open("%sURL.txt" % (txt_name),'w')
for i,line in enumerate(fIN):
    fOUT.write(site_redirector[site]+line)
fOUT.close()
fIN.close()