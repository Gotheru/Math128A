function [root,info] = modifiedzeroin(func, Int, params)

    a = Int.a;
    b = Int.b;
    fa = func(a);
    fb = func(b);
    assert(a <= b, "Invalid input: a > b")
    assert(fa*fb <= 0, "Invalid input: f(a) * f(b) > 0");
    x0 = Int.a;
    f0 = fa;
    x1 = Int.b;
    f1 = fb;
    x2 = (Int.a + Int.b) / 2;
    f2 = func(x2);

    %% We store f0, f1, f2 so that we don't make as many function calls.
    %% fa and fb for the same reason.
    
    % Define lambda function iqi
    iqi = @(x0, f0, x1, f1, x2, f2) x0 * f1 * f2 / (f0 - f1) / (f0 - f2) + x1 * f0 * f2 / (f1 - f0) / (f1 - f2) + x2 * f0 * f1 / (f2 - f0) / (f2 - f1);
    
    function_calls = 3;

    info.flag = 1;

    consecutive_bad_iterations = 0; % Counter of consecutive bad iterations.
    
    for it = 1:params.maxit
        x3 = iqi(x0, f0, x1, f1, x2, f2);
        f3 = func(x3);
        function_calls = function_calls + 1;

        if 2 * abs(f3) > abs(f2)
            consecutive_bad_iterations = consecutive_bad_iterations + 1;
        else
            consecutive_bad_iterations = 0;
        end

        if x3 < a || x3 > b || consecutive_bad_iterations >= 4
            % do bisection steps on [a, b]
            c = (a + b) / 2;
            fc = func(c);
            function_calls = function_calls + 1;
            if fa * fc <= 0
                b = c;
                fb = fc;
            else
                a = c;
                fa = fc;
            end
            x0 = a;
            f0 = fa;
            x1 = b;
            f1 = fb;
            x2 = (a + b) / 2;
            f2 = func(a + b) / 2;
            function_calls = function_calls + 1;
            % One less bad iteration (adjustment)
            consecutive_bad_iterations = max(0, consecutive_bad_iterations - 1);
        else
            % Iterate
            x0 = x1;
            f0 = f1;
            x1 = x2;
            f1 = f2;
            x2 = x3;
            f2 = f3;
            % Adjust a and b: make a smaller interval so we make less
            % calls.
            if f2 * fa <= 0
                b = x2;
                fb = f2;
            else
                a = x2;
                fa = f2;
            end
        end
        
        if abs(f2) < params.func_tol || b - a < params.root_tol
            info.flag = 0;
            info.iterations = it;
            info.function_calls = function_calls;
            break;
        end
        
    end
    % Best approximation (or the actual value if we break)
    root = x2;
    info.function_calls = function_calls;

end