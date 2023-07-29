#!/bin/bash
#2020-02-17
input=$1

#define run settings
exe_path="~/rosetta/main/source/bin/rosetta_scripts.macosclangrelease"
inpath=/sandbox/relax/${input}
outpath=output/${input}

score_function="beta"
    if [[ ${score_function} == "beta" ]]; then beta="-beta"; fi

if [ ! -e output/ ]; then mkdir output/; fi
if [ ! -e output/${input} ]; then mkdir output/${input}; fi

outpath=output/${input}

#run Rosetta
for pdb in $(ls ${inpath}/*.pdb);do
	in_name=`echo ${pdb##*/} |cut -d'.' -f1`
	${exe_path} \
		${beta} \
        	-dunbrack_prob_buried 0.8 -dunbrack_prob_nonburied 0.8 -dunbrack_prob_buried_semi 0.8 -dunbrack_prob_nonburied_semi 0.8 \
		-out::file::pdb_comments \
		-parser:protocol xml/filter_scan_int.xml \
		-s      ${pdb} \
		-native ${pdb} \
		-nstruct 1 \
		-overwrite 1 \
		-unmute all \
		-out:chtimestamp 1 \
		-out:suffix "" \
		-out:level 300 \
		-out::path::all ${outpath}/ \
		-failed_job_exception False \
		-mute core.select.residue_selector.SecondaryStructureSelector \
		> ${outpath}/${in_name}.log
done
