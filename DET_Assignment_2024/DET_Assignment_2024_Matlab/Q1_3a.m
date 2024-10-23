function data = generate_exponential_data(num_datapoints, lambda)
    % Generates i.i.d. samples from an exponential distribution using the
    % inverse transform sampling method.
    % Inputs:
    % - num_datapoints: Number of datapoints to generate
    % - lambda: The true rate parameter of the exponential distribution
    % Output:
    % - data: A vector of generated data points

    % Generate uniform random numbers between 0 and 1
    U = rand(num_datapoints, 1);
    
    % Transform uniform random numbers to exponential random numbers
    data = -log(1 - U) / lambda;
end

function crlb = compute_crlb(data,lambda)
    % Computes the Cramer Rao Lower Bound
    % Inputs:
    % - data: A vector of observed data points from an exponential distribution
    % - lambda : The true value of lambda
    % Output:
    % - crlb : The Cramer Rao Lower Bound
    
    % Compute the CRLB given from λ^2 / Ν
    N = length(data);
    crlb = (lambda)^2 / N;
end

function lambda_mvue = compute_mvue(data)
    % Computes an MVUE that attains the CRLB for lambda
    % Inputs:
    % - data: A vector of observed data points from an exponential distribution
    % - lambda : The true value of lambda
    % Output:
    % - lambda_mvue: The estimated rate parameter (lambda)
   
    
    % The MVUE for lambda in an exponential distribution is
    % n / sum(X)

    n = length(data);
    lambda_mvue = n / sum(data);
    
end

function mvue_results = perform_mvue_experiments(num_experiments, num_datapoints, lambda)
    % Performs multiple experiments and computes an MVUE for each experiment
    % Inputs:
    % - num_experiments: Number of experiments to perform
    % - num_datapoints: Number of datapoints per experiment
    % - lambda: The true rate parameter of the exponential distribution
    % Output:
    % - mvue_results: A vector of MVUE estimates for each experiment

    mvue_results = zeros(num_experiments, 1);
    
    for i = 1:num_experiments
        data = generate_exponential_data(num_datapoints, lambda);
        mvue_results(i) = compute_mvue(data);
    end
end


% Main script to run the entire process
lambda = 2;
number_datapoints_1 = 20;
number_datapoints_2 = 50;
number_datapoints_3 = 100;
number_experiments = 5000;

% Generate a random dataset for the CRLB calculations
data_1 = generate_exponential_data(number_datapoints_1,lambda);
data_2 = generate_exponential_data(number_datapoints_2,lambda);
data_3 = generate_exponential_data(number_datapoints_3,lambda);

% Calculate and display the CRLBs
crlb_1 = compute_crlb(data_1,lambda);
crlb_2 = compute_crlb(data_2,lambda);
crlb_3 = compute_crlb(data_3,lambda);
disp(['Ν = 20 : CRLB at lambda = ', num2str(lambda), ': ', num2str(crlb_1)]);
disp(['Ν = 50 : CRLB at lambda = ', num2str(lambda), ': ', num2str(crlb_2)]);
disp(['Ν = 100 : CRLB at lambda = ', num2str(lambda), ': ', num2str(crlb_3)]);

% Compute the Unbiased Estimators using the function
mvue_1 = perform_mvue_experiments(number_experiments,number_datapoints_1,lambda);
mvue_2 = perform_mvue_experiments(number_experiments,number_datapoints_2,lambda);
mvue_3 = perform_mvue_experiments(number_experiments,number_datapoints_3,lambda);

% Plot histogram of MVUE estimates for N = 20
figure;
histogram(mvue_1, 50);
title('Histogram of Unbiased Estimates(N=20)');
xlabel('Estimated \lambda');
ylabel('Frequency');

% Plot histogram of MVUE estimates for N = 50
figure;
histogram(mvue_2, 50);
title('Histogram of Unbiased Estimates(N=50)');
xlabel('Estimated \lambda');
ylabel('Frequency');

% Plot histogram of MVUE estimates for N = 100
figure;
histogram(mvue_3, 50);
title('Histogram of Unbiased Estimates(N=100)');
xlabel('Estimated \lambda');
ylabel('Frequency');