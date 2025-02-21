addpath('../generic');

%%

ndatasets = 8;

dataset = [];
choice = [];
RT = [];
correct = [];
coh = [];
coh_extra = [];
monkey = [];
pulse_size = [];
pulse_on = [];
pulse_dur = [];
for idataset = 1:ndatasets

    disp(num2str(idataset));

    d = get_data(nan,[],idataset,'positive_is_leftward',1);

    ntr = length(d.choice);

    dataset = [dataset; idataset*ones(ntr,1)];
    choice  = [choice; d.choice];
    RT      = [RT; d.RT];
    correct = [correct; d.correct];
    coh     = [coh; d.coh];
    coh_extra = [coh_extra; d.coh_extra];
    monkey = [monkey; d.monkey];

    pulse_on = [pulse_on; d.pulse_on];
    pulse_dur = [pulse_dur; d.pulse_off - d.pulse_on];
    pulse_size = [pulse_size; d.pulse_size/1000]; 

end

[uni_monkey,~,monkey_id] = unique(monkey);

save('data_for_DTB','dataset','choice','RT','correct','coh','coh_extra','monkey_id','uni_monkey','pulse_on','pulse_dur','pulse_size'); % no H


