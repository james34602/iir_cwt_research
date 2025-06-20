function eq = fit2eq(cf)
eq = formula(cf); %Formula of fitted equation

parameters = coeffnames(cf); %All the parameter names

values = coeffvalues(cf); %All the parameter values

for idx = 1:numel(parameters)

    param = parameters{idx};

    l = length(param);

    loc = regexp(eq, param); %Location of the parameter within the string

    while ~isempty(loc)

        %Substitute parameter value

        eq = [eq(1:loc-1) num2str(values(idx), 14) eq(loc+l:end)];

        loc = regexp(eq, param);

    end

end
end