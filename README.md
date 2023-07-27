# Documentation

#https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=37881&filename=xsec_tools_analysis_retreat_may2022.pdf&version=3

#https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=37881&filename=documentation_tech_note_partial_draft.pdf&version=3

#https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=35518&filename=internal_note_v04.pdf&version=9 


cd /uboone/app/users/USERNAME
git clone https://github.com/afropapp13/xsec_analyzer.git ./
cd xsec_analyzer
source setup.sh

# To define the event category (backgrounds et al)
emacs -nw EventCategory.hh

# FilePropertiesManager.hh defines the allowed file types, etc
emacs -nw FilePropertiesManager.hh

make

# Generate the configuration
root -b make_config_all.C

# Pro-processing
./reprocess.sh /uboone/data/users/apapadop/PeLEE_PostProcessing/   files_to_process.txt

# Dumping the post-processed file paths in txt file
ls -d /uboone/data/users/apapadop/PeLEE_PostProcessing/*.root | tee post_processed_file_list.txt

# Systematics manipulation
./univmake post_processed_file_list.txt myconfig_all.txt /uboone/data/users/apapadop/PeLEE_PostProcessing/xsec_analysis.root file_properties.txt

