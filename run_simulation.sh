# Spatial dependence
Rscript src/3_model_usage_and_simulation/1_estimate_spatial_dependence_structure/2_seasonal.R '1_parametric' 'glm'
Rscript src/3_model_usage_and_simulation/1_estimate_spatial_dependence_structure/2_seasonal.R '1_parametric' 'bamlss'

# Simulation
Rscript src/3_model_usage_and_simulation/2_simulate/2_seasonal_joint_latent_scheme.R '1_parametric' 'glm'
Rscript src/3_model_usage_and_simulation/2_simulate/2_seasonal_joint_latent_scheme.R '1_parametric' 'bamlss'

