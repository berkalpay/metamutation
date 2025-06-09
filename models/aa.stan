data {
  int<lower=1> P;
  int<lower=1> N[P];
  int<lower=1> maxN;
  int<lower=0, upper=20> w[P, maxN];
  int<lower=0, upper=20> m[P, maxN];
}
parameters {
  matrix[20, 20] mu;
}
model {
  matrix[maxN, P] score;
  for (wi in 1:20) {
    for (mi in 1:20) {
      mu[wi, mi] ~ normal(0, 1);
    }
  }
  for (i in 1:P) {
    for (j in 1:N[i]) {
      score[j, i] = 1/(1+exp(-mu[w[i, j], m[i, j]]));
    }
    vector[N[i]] scores = score[1:N[i], i];
    target += sum(log(scores) - log(reverse(cumulative_sum(reverse(scores)))));
  }
}
