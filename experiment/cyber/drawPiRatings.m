
Alpha = runner.combiner{2}.Alpha;

Nu = runner.combiner{2}.Nu;%runner.combiner{2}.nu + sum(results, 2)';

[ratings, classRatings, igRatings] = infgain.piRating(Alpha, Nu);
figure;
bar(sort(ratings));
xlabel('Features');
ylabel('Information value');

figure;
bar(sort(igRatings));
xlabel('Features');
ylabel('Information value');