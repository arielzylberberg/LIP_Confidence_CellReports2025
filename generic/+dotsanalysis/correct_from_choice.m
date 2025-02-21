function correct = correct_from_choice(choice, coh)

correct = nan(size(choice));
I = coh>0;
correct(I) = choice(I);
I = coh<0;
correct(I) = 1-choice(I);

I = coh==0;
correct(I) = 0.5; % clamp it

end