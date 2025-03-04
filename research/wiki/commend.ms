snakemake all --configfile=config/config.public.yaml

snakemake ariadne_all --configfile=config/config.public.yaml

snakemake write_statistics --configfile=config/config.public.yaml
snakemake make_summary --configfile=config/config.public.yaml
snakemake -call plot_statistics_comparison --configfile=config/config.public.yaml


