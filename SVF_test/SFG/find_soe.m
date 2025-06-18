function out = find_soe(tfq,TFqs,numonly)

[NUM,DEN]=tfdata(tf(tfq),'v');
[NUMs,DENs]=numden((TFqs));

NUMs = collect(NUMs,'z');
DENs = collect(DENs,'z');

denfirst = 1;

mask = find(abs(NUM)>0);

if not(denfirst)
    [out1,MUL] = extract_sym_tf(NUMs,NUM(mask(1)));
    out2 = extract_sym_tf(MUL*DENs,DEN(1));
else
    [out2,MUL] = extract_sym_tf(DENs,DEN(1));
    out1 = extract_sym_tf(MUL*NUMs,NUM(mask(1)));
end

if nargin == 2
    out = [out1 ;out2];
    NNN = [NUM.';DEN.'];
    ttt = zeros(size(NNN));
elseif nargin == 3
    if numonly == 1
        out = out1;
        NNN = NUM.';
        ttt = zeros(size(NNN));
    end
end

for ind = 1:length(out)
    
    temp = symvar(out(ind));
    ttt(ind) = length(temp);
end

out=out(ttt>0);
NNN=NNN(ttt>0);

L=length(NNN);
for ind = 1:L
    symnnn=sym(NNN(ind));
    out(ind) = str2sym([char(out(ind)),'-(',char(symnnn),')']);
end

end

function varargout = extract_sym_tf(ND,mul)

if nargin == 1
    mul=1;
end

K = subs(ND,'z',0);
ND = ND-K;
temp = char(collect(ND,'z'));

mask = strfind(temp,'z');
pot = zeros(size(mask));

L=length(temp);

for ind = length(mask):-1:1

    if mask(ind)+1 < L
        if strfind(temp(mask(ind)+1),'^')==1
            pot(ind) = eval(temp(mask(ind)+2));
        else
            temp = [temp(1:mask(ind)),'^1',temp((mask(ind)+1):end)];
            pot(ind) = 1;
        end
    else
        temp = [temp,'^1'];
        pot(ind) = 1;
    end
end

mask = strfind(temp,'z');
%A = sym('A',[length(pot),1]);
A = sym(zeros(length(pot),1));
for ind = 1:length(pot)
    myeq = temp;
    pmask = mask(not(pot==ind));
    for indm = 1:length(pmask)
        myeq(pmask(indm):(pmask(indm)+2))='0.0';
    end
    A(ind) = simplify(str2sym(['(',myeq,')/',myeq(mask(pot==ind):mask(pot==ind)+2)]));
end


A = flipud(A);
A=[A;K];

Am=1/A(1);

A=mul*A*Am;

varargout{1}=A;

if nargout == 2
    varargout{2}=Am;
end

end