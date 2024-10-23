function data = generate_exponential_data(num_datapoints, true_lambda)
    % Generates i.i.d. samples from an exponential distribution using the
    % inverse transform sampling method.
    % Inputs:
    % - num_datapoints: Number of datapoints to generate
    % - true_lambda: The true rate parameter of the exponential distribution
    % Output:
    % - data: A vector of generated data points

    % Generate uniform random numbers between 0 and 1
    U = rand(num_datapoints, 1);
    
    % Transform uniform random numbers to exponential random numbers
    data = -log(1 - U) / true_lambda;
end

function lambda_mle = compute_mle(data)
    % Computes the Maximum Likelihood Estimator (MLE) for lambda
    % Inputs:
    % - data: A vector of observed data points from an exponential distribution
    % Output:
    % - lambda_mle: The estimated rate parameter (lambda)
    
    % The MLE for lambda in an exponential distribution is 1 / mean(data)
    lambda_mle = 1 / mean(data);
end

function mle_results = perform_mle_experiments(num_experiments, num_datapoints, lambda)
    % Performs multiple experiments and computes MLE for each experiment
    % Inputs:
    % - num_experiments: Number of experiments to perform
    % - num_datapoints: Number of datapoints per experiment
    % - lambda: The true rate parameter of the exponential distribution
    % Output:
    % - mle_results: A vector of MLE estimates for each experiment

    mle_results = zeros(num_experiments, 1);
    
    for i = 1:num_experiments
        data = generate_exponential_data(num_datapoints, lambda);
        mle_results(i) = compute_mle(data);
    end
end


% Generate example datasets
number_datapoints_1 = 20;
number_datapoints_2 = 50;
number_datapoints_3 = 100;
lambda = 2;
number_experiments = 5000;

% Compute MLE estimator for all cases
lambda_mle_1 = perform_mle_experiments(number_experiments,number_datapoints_1,lambda);
lambda_mle_2 = perform_mle_experiments(number_experiments,number_datapoints_2,lambda);
lambda_mle_3 = perform_mle_experiments(number_experiments,number_datapoints_3,lambda);

% Plot histogram of MVUE estimates for N = 20
figure;
histogram(lambda_mle_1, 50);
title('Histogram of MLE Estimates(N = 20)');
xlabel('Estimated \lambda');
ylabel('Frequency');

% Plot histogram of MVUE estimates for N = 50
figure;
histogram(lambda_mle_2, 50);
title('Histogram of MLE Estimates(N = 50)');
xlabel('Estimated \lambda');
ylabel('Frequency');

% Plot histogram of MVUE estimates for N = 100
figure;
histogram(lambda_mle_3, 50);
title('Histogram of MLE Estimates(N = 100)');
xlabel('Estimated \lambda');
ylabel('Frequency');