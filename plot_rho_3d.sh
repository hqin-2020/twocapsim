#! /bin/bash

# Deltaarray=(100 500 1000)
# Deltaarray=(50 100 150 300)
Deltaarray=(100 150)
fractionarray=(0.0)

actiontime=1

python_name="plots_rho_moregrid.py"

rhoarray=(0.7 0.8 0.9 1.00001 1.1 1.2 1.3 1.4 1.5)
# rhoarray=(1.00001)
gammaarray=(8.0)
# gammaarray=(1.01 1.1 1.25 1.5 2.0 4.0 8.0)
# gammaarray=(4.1 4.3 4.5 5.0 6.0)
# gammaarray=(4.6 4.7 4.8 4.9)
# gammaarray=(1.1 2.0 4.0 5.0 8.0)

for Delta in ${Deltaarray[@]}; do
    for fraction in "${fractionarray[@]}"; do
        for rho in "${rhoarray[@]}"; do
            for gamma in "${gammaarray[@]}"; do
                    count=0

                    action_name="TwoCapital_julia_rhoeq_required_test_more_grid4"
                    action_name="TwoCapital_julia_rhoeq_more_grids_gamma_56"
                    action_name="TwoCapital_julia_rhoeq_standard_grids0"
                    # action_name="TwoCapital_julia_rhoeq_standard_grids_as"

                    dataname="${action_name}_${Delta}_frac_${fraction}"

                    mkdir -p ./job-outs/${action_name}/Delta_${Delta}_frac_${fraction}/

                    if [ -f ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.sh ]; then
                        rm ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.sh
                    fi

                    mkdir -p ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/

                    touch ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.sh

                    tee -a ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=p_${Delta}_${fraction}
#SBATCH --output=./job-outs/$job_name/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.out
#SBATCH --error=./job-outs/$job_name/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.err
#SBATCH --time=0-12:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G

module load python/anaconda-2020.11
python3 /project/lhansen/twocapsim/$python_name  --Delta ${Delta} --fraction ${fraction} --gamma ${gamma} --rho ${rho} --dataname ${dataname}
EOF
                count=$(($count + 1))
                sbatch ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}_plot.sh
            done
        done
    done
done