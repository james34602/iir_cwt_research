classdef ADNode2 < handle
    properties
        value % function value at this node
        grad % gradient accumulator
        func % callback function to update gradient of the parent nodes
        root % input node that holds the tape
        tape % sequence of evaluation steps
        str0
        str1
        str2
    end
    methods
        function y = ADNode2(str0, str1, str2, x, root, func)
            if nargin > 4
                y.func = func;
                y.root = root;
                root.tape{end+1} = y;
            else
                y.root = y;
                y.tape = {};
            end
            y.str0 = str0;
            y.str1 = str1;
            y.str2 = str2;
            y.value = x;
        end
        function changeName(obj, name)
            obj.str0 = string(name);
        end
        function y = sum(x, dim, flag)
            switch nargin
              case 3
                y = ADNode2('','','', sum(x.value, dim, flag), x.root, @(y) x.add(y.grad .* ones(size(x.value))));
              case 2
                y = ADNode2('','','', sum(x.value, dim), x.root, @(y) x.add(y.grad .* ones(size(x.value))));
              otherwise
                y = ADNode2('','','', sum(x.value), x.root, @(y) x.add(y.grad .* ones(size(x.value))));
            end
        end
        function y = reshape(x, shp)
            y = ADNode2('','','', reshape(x.value, shp), x.root, @(y) x.add(reshape(y.grad, size(x.value))));
        end
        function y = cos(x)
            y = ADNode2('','','', cos(x.value), x.root, @(y) x.add(bsxfun(@times, -y.grad, sin((x.value)))));
        end
        function dy = backprop(x)
            x.grad = 1;
            % progress, k = 91
            for k = length(x.root.tape):-1:1
                if k == 91
                    disp('');
                end
                x.root.tape{k}.func(x.root.tape{k});
                x.root.tape(k) = [];
            end
            dy = x.root.grad;
        end
        function [varargout] = subsref(x, s)
            switch s(1).type
              case '()'
                varargout{1} = ADNode2(x.str0, x.str1, x.str2, x.value(s.subs{:}), x.root, @(y) x.subs_add(s.subs, y));
              otherwise
                  [varargout{1:nargout}] = builtin('subsref', x, s);
            end
        end
        function y = subsasgn(x, s, varargin)
            inputOrigSize = size(x.value);
            x.value(s.subs{:}) = varargin{1}.value;
            t = ADNode2('','','', x.value(s.subs{:}), x.root, @(y) varargin{1}.subs_move(s.subs, x, inputOrigSize));
            y = x;
        end
        function y = plus(x1, x2, str0, str1, str2)
            disp(str0 + '_value' + ' = ' + str1 + '_value' + ' + ' + str2 + '_value')
            y = ADNode2(str0, str1, str2, bsxfun(@plus, x1.value, x2.value), x1.root, @(y) y.plus_backprop(x1, x2, str0, str1, str2));
        end
        function y = minus(x1, x2, str0, str1, str2)
            disp(str0 + '_value' + ' = ' + str1 + '_value' + ' - ' + str2 + '_value')
            y = ADNode2(str0, str1, str2, bsxfun(@minus, x1.value, x2.value), x1.root, @(y) y.minus_backprop(x1, x2, str0, str1, str2));
        end
        function aa(x1, x2, y)
            x1.add(y.grad * x2');
        end
        function y = mtimes(x1, x2)
            y = ADNode2('','','', x1.value * x2, x1.root, @(y) aa(x1, x2, y));
        end
        function y = times(x1, x2, str0, str1, str2)
            disp(str0 + '_value' + ' = ' + str1 + '_value' + ' * ' + str2 + '_value')
            y = ADNode2(str0, str1, str2, x1.value * x2.value, x1.root, @(y) y.times_backprop(x1, x2, str0, str1, str2));
        end
        function y = power(x1, x2)
            y = ADNode2('','','',x1.value .^ x2, x1.root, @(y) x1.add(bsxfun(@times, y.grad, x1.value * 2)));
        end
        function y = length(adn)
            y = length(adn.value);
        end
        function y = size(adn, dim)
            if nargin < 2
                y = size(adn.value);
            else
                y = size(adn.value, dim);
            end
        end
    end
    methods (Access = private)
        function add(x, grad)
            %% accumulate the gradient, take sum of dimensions if needed
            if isempty(x.grad)
                if size(x.value) == [1, 1]
                    x.grad = sum(sum(grad));
                elseif size(x.value, 1) == 1
                    x.grad = sum(grad, 1);
                elseif size(x.value, 2) == 1
                    x.grad = sum(grad, 2);
                else
                    x.grad = grad;
                end
            else
                if size(x.grad) == [1, 1]
                    x.grad = x.grad + sum(sum(grad));
                elseif size(x.grad, 1) == 1
                    x.grad = x.grad + sum(grad, 1);
                elseif size(x.grad, 2) == 1
                    x.grad = x.grad + sum(grad, 2);
                else
                    x.grad = bsxfun(@plus, x.grad, grad);
                end
            end
        end
        function transp(x, grad)
            x.grad = grad';
        end

        function subs_add(x, subs, y)
            %% accumulate the gradient with subscripts
            gradt = y.grad;
            if isempty(x.grad)
                x.grad = zeros(size(x.value));
            end
            old = x.grad(subs{:});
            if size(old, 1) == 1 && size(old, 2) == 1
                x.grad(subs{:}) = old + sum(sum(gradt));
            elseif size(old, 1) == 1
                x.grad(subs{:}) = old + sum(gradt, 1);
            elseif size(old, 2) == 1
                x.grad(subs{:}) = old + sum(gradt(1 : min(size(old, 1), size(gradt, 1)), :), 2);
            else
                x.grad(subs{:}) = old + gradt;
            end
        end

        function subs_match(x, subs, y)
            %% accumulate the gradient with subscripts
            if isempty(x.grad)
                x.grad = zeros(size(x.value));
            end
            if size(x.grad) == [1, 1]
                x.grad = x.grad + sum(y.grad(subs{:}));
            else
                x.grad(subs{:}) = x.grad(subs{:}) + y.grad(subs{:});
            end
        end

        function subs_clear(x, subs)
            %% clear the gradient with subscripts
            if isempty(x.grad)
                x.grad = zeros(size(x.value));
            end
            x.grad(subs{:}) = 0;
        end

        function subs_move(x, subs, y, inputOrigSize)
            %% accumulate the gradient with subscripts
            if size(y.grad) == [1,1]
                y.grad = repmat(y.grad, size(y.value));
            end
            gradt = y.grad(subs{:});
            y.grad(subs{:}) = 0;
            if any(size(y.grad) ~= inputOrigSize) % Remove excessive gradient data
                y.grad(subs{:}) = [];
            end
            if isempty(x.grad)
                x.grad = zeros(size(x.value));
            end
            old = x.grad;
            if size(old, 1) == 1 && size(old, 2) == 1
                x.grad = old + sum(sum(gradt));
            elseif size(old, 1) == 1
                x.grad = old + sum(gradt, 1);
            elseif size(old, 2) == 1
                x.grad = old + sum(gradt, 2);
            else
                x.grad = old + gradt;
            end
        end
        function plus_backprop(y, x1, x2, str0, str1, str2)
            disp(x1.str0 + '_grad = ' + x1.str0 + '_grad + ' + y.str0 + '_grad');
            disp(x2.str0 + '_grad = ' + x2.str0 + '_grad + ' + y.str0 + '_grad');
            if contains(x2.str0, "wd2")
                disp('');
            end
            x1.add(y.grad);
            x2.add(y.grad);
        end
        function minus_backprop(y, x1, x2, str0, str1, str2)
            disp(x1.str0 + '_grad = ' + x1.str0 + '_grad + ' + y.str0 + '_grad');
            disp(x2.str0 + '_grad = ' + x2.str0 + '_grad - ' + y.str0 + '_grad');
            if contains(x2.str0, "z1_A")
                disp('');
            end
            x1.add(y.grad);
            x2.add(-y.grad);
        end
        function times_backprop(y, x1, x2, str0, str1, str2)
            disp(x1.str0 + '_grad = ' + x1.str0 + '_grad + ' + y.str0 + '_grad * ' + x2.str0 + '_value');
            disp(x2.str0 + '_grad = ' + x2.str0 + '_grad + ' + y.str0 + '_grad * ' + x1.str0 + '_value');
%             if (x2.str0 == "z2_A(4)")
%                 disp('');
%             end
            if contains(x1.str0, "c1")
                disp('');
            end
            if contains(x2.str0, "z1_A")
                disp('');
            end
            x1.add(bsxfun(@times, y.grad, x2.value));
            x2.add(bsxfun(@times, y.grad, x1.value));
        end
        function mtimes_backprop(y, x1, x2)
            x1.add(y.grad * x2.value');
            x2.add(x1.value' * y.grad);
        end
    end
end