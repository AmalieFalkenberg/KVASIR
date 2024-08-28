for d_AD in $(seq 1 0.15 4); do 
            mcrun -c --mpi=6 KVASIR1p0.instr -d FWHM_E_d_AD_${d_AD} -n 1000000 -N 41 Ef=1.83 d_SA=1,4 d_AD=${d_AD} dtA=0;
done

for d_AD in $(seq 1 0.15 4); do 
            mcrun -c --mpi=6 KVASIR1p0.instr -d FWHM_E19_d_AD_${d_AD} -n 1000000 -N 41 Ef=1.9 d_SA=1,4 d_AD=${d_AD} dtA=0;
done

for d_AD in $(seq 1 0.15 4); do 
            mcrun -c --mpi=6 KVASIR1p0.instr -d FWHM_E2_d_AD_${d_AD} -n 1000000 -N 41 Ef=2 d_SA=1,4 d_AD=${d_AD} dtA=0;
done