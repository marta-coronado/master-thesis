#! bin/sh

ls -tr ../raw_data/LIST90 > LIST90
for file in `less LIST90`
	do
		perl input.pl -file ../raw_data/LIST90/${file} #it generates the inputs > ../dfe-alpha-release-2.13/sfs.txt and ../dfe-alpha-release-2.13/divergence.txt
		sed  's/:/ /g' ./INPUT/sfs.tmp > ./INPUT/sfs.txt
		rm ./INPUT/sfs.tmp
		./est_dfe -c ./CONFIG/config-file-dfe-0.txt #to estimate demography with neutral sites
		./est_dfe -c ./CONFIG/config-file-dfe-1.txt #to estimate the DFE of selective sites taking into account the demography infered with the neutral sites
		./est_alpha_omega -c ./CONFIG/config-file-alpha_omega.txt #to estimate alpha taking into account all the previous estimations (demography + DFE)
		# ./est_alpha_omega -c ./CONFIG/example-config-file-for-est_alpha_omega.txt #>> ./OUTPUT/results_alpha_sel/history
		
		mv ./OUTPUT/results_alpha_sel/est_alpha_omega.out  ../results_order2cat_ex/${file}.alpha.txt
		echo "${file}" >> ../results_order2cat_ex/RESULTS_ALPHA_OMEGA
		cat ../results_order2cat_ex/${file}.alpha.txt >> ../results_order2cat_ex/RESULTS_ALPHA_OMEGA

		mv ./OUTPUT/results_dir_sel/est_dfe.out  ../results_order2cat_ex/${file}.dfe.txt
		echo "${file}" >> ../results_order2cat_ex/RESULTS_DFE
		cat ../results_order2cat_ex/${file}.dfe.txt >> ../results_order2cat_ex/RESULTS_DFE
		
		rm ../results_order2cat_ex/${file}.alpha.txt ../results_order2cat_ex/${file}.dfe.txt
	done	


exit 0

<<COMMENT

ls -tr ../raw_data/LIST > LIST

for file in `less LIST`
	do
		perl input.pl -file ../raw_data/LIST/${file} #it generates the inputs > ../dfe-alpha-release-2.13/sfs.txt and ../dfe-alpha-release-2.13/divergence.txt
		sed  's/:/ /g' ./INPUT/sfs.tmp > ./INPUT/sfs.txt
		rm ./INPUT/sfs.tmp
		./est_dfe -c ./CONFIG/config-file-dfe-0.txt #to estimate demography with neutral sites
		./est_dfe -c ./CONFIG/config-file-dfe-1.txt #to estimate the DFE of selective sites taking into account the demography infered with the neutral sites
		./est_alpha_omega -c ./CONFIG/config-file-alpha_omega.txt #to estimate alpha taking into account all the previous estimations (demography + DFE)
		# ./est_alpha_omega -c ./CONFIG/example-config-file-for-est_alpha_omega.txt #>> ./OUTPUT/results_alpha_sel/history
		
		mv ./OUTPUT/results_alpha_sel/est_alpha_omega.out  ../results_dfe/${file}.alpha.txt
		echo "${file}" >> ../results_dfe/RESULTS_ALPHA_OMEGA
		cat ../results_dfe/${file}.alpha.txt >> ../results_dfe/RESULTS_ALPHA_OMEGA

		mv ./OUTPUT/results_dir_sel/est_dfe.out  ../results_dfe/${file}.dfe.txt
		echo "${file}" >> ../results_dfe/RESULTS_DFE
		cat ../results_dfe/${file}.dfe.txt >> ../results_dfe/RESULTS_DFE
		
		rm ../results_dfe/${file}.alpha.txt ../results_dfe/${file}.dfe.txt
	done	


exit 0
COMMENT