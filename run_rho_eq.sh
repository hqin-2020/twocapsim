#! /bin/bash

# Deltaarray=(5 10 20 30 50 100 150 500 1000 5000 10000)
# fractionarray=(0.1 0.05 0.01 0.005 0.001 0.0005 0.0001)
# Deltaarray=(5 10 20 50 100 150 300)
Deltaarray=(5 50 100 150 300)
fractionarray=(0.0)

actiontime=1

julia_name="newsets_twocapitals_rhoas0.jl"

# rhoarray=(0.7 0.8 0.9 1.1 1.2 1.3 1.4 1.5)
rhoarray=(1.00001)
# gammaarray=(1.1 2.0 4.0 5.0 8.0)
# gammaarray=(8.0)
# gammaarray=(1.01 1.05 3.0 5.0)
# gammaarray=(1.00001 2.0 4.0 8.0)
gammaarray=(1.01 1.1 1.5 2.0 4.0 8.0)
# gammaarray=(4.1 4.3 4.5 5.0 6.0)
# gammaarray=(4.6 4.7 4.8 4.9)

for Delta in ${Deltaarray[@]}; do
    for fraction in "${fractionarray[@]}"; do
        for rho in "${rhoarray[@]}"; do
            for gamma in "${gammaarray[@]}"; do
                    count=0

                    # action_name="TwoCapital_julia_rhoeq_more_test"
                    # action_name="TwoCapital_julia_rhoeq_more_grids_gamma_56"
                    action_name="TwoCapital_julia_rhoeq_standard_grids0"
                    action_name="TwoCapital_julia_rhoas_standard_grids0"

                    dataname="${action_name}_${Delta}_frac_${fraction}"

                    mkdir -p ./job-outs/${action_name}/Delta_${Delta}_frac_${fraction}/

                    if [ -f ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.sh ]; then
                        rm ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.sh
                    fi

                    mkdir -p ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/

                    touch ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.sh

                    tee -a ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=${Delta}_${rho}
#SBATCH --output=./job-outs/$job_name/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.out
#SBATCH --error=./job-outs/$job_name/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.err
#SBATCH --time=0-12:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module load julia/1.7.3
srun julia /project/lhansen/twocapsim/$julia_name  --Delta ${Delta} --fraction ${fraction} --gamma ${gamma} --rho ${rho} --dataname ${dataname}
EOF
                count=$(($count + 1))
                sbatch ./bash/${action_name}/Delta_${Delta}_frac_${fraction}/rho_${rho}_gamma_${gamma}.sh
            done
        done
    done
done