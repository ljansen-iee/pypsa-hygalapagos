snakemake solve_network --cluster-config config.cluster.yaml \
--cluster "sbatch -p {cluster.partition} -t {cluster.walltime} -c {cluster.cpus_per_task} --mem {cluster.mem_mb} -x {cluster.exclude}" \
--jobs 1000 --latency-wait 60 --keep-going \
--rerun-trigger code \
# --forceall \
# --rerun-incomplete \
# -n \


# squeue --sort=t,p,-S --format="%.7i %.9P %.8u %.8j %.4Q %.7T %.19V %.19S %.11M %.11l %.19e %.2D %.3C %R" -p progress

# sinfo -o "%n %e %m %a %c %C" | sort -k5,5nr