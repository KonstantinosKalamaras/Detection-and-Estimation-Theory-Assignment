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
number_datapoints = 10;
number_experiments = 5000;

% Generate a random dataset for the CRLB calculation
data = generate_exponential_data(number_datapoints,lambda);

% Calculate and display the CRLB
crlb = compute_crlb(data,lambda);
disp(['CRLB at lambda = ', num2str(lambda), ': ', num2str(crlb)]);

% Compute the MVUE using the function
mvue = perform_mvue_experiments(number_experiments,number_datapoints,lambda);


% Plot histogram of MVUE estimates
figure;
histogram(mvue, 50);
title('Histogram of MVUE Estimates');
xlabel('Estimated \lambda');
ylabel('Frequency');