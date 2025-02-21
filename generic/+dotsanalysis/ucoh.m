function ucoh = ucoh()
ucoh = [0,0.032*2.^[0:4]];
ucoh = sort(unique([-ucoh,ucoh]));