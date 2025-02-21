load('./output_data/out_boot_shuffle');


%%

mean(fit_alpha_boot(1,:)>fit_alpha_boot(2,:))
mean(fit_alpha_boot(2,:)>fit_alpha_boot(3,:))

