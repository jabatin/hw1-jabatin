% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
%x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
 x = 3; y = '5'; %mixed

%your code goes here
a = ischar(x);
b = ischar(y);
if (a == 1) 
    c = str2num(x);
    x = c;
end
if (b == 1) 
    d = str2num(y);
    y = d;
end
%output your answer
z = x + y;
z

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length
Index = 'AGCT';
rand_seq = '';
count = 0
while count < 500
    Vari = Index(randi(numel(Index)));
    rand_seq = strcat(rand_seq, Vari);
    count = numel(rand_seq);
end
rand_seq
%%
%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

SeqSt = strfind(rand_seq, 'ATG');
SeqEd1 = strfind(rand_seq, 'TAA');
SeqEd2 = strfind(rand_seq, 'TGA');
SeqEd3 = strfind(rand_seq, 'TAG');
SeqEdAll = [SeqEd1, SeqEd2, SeqEd3];
SeqEd = sort(SeqEdAll);
x = numel(SeqSt);
b = numel(SeqEd);
y = 1;
SeqNew = [];
SeqFinal = [];
while y<(x+1)
    a = 1;
    while a<(b+1)
    Sub = SeqSt(y)- SeqEd(a);
    SeqNew = [SeqNew, Sub];
    a = a + 1;
    end
    Seq2 = SeqNew(SeqNew>=0);
    Seq2 = sort(Seq2);
    if numel(Seq2) < 1
    SeqFinal = SeqFinal
    else
    SeqFinal = [SeqFinal, Seq2(1)]
    end
    y = y + 1;
end
SeqFinal = fliplr(sort(SeqFinal));
if numel(SeqFinal) < 1
LongORF = 'undefined'
else
LongORF = SeqFinal(1);
end
%%
%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.
counter = 1;
ProbArray = [];
while counter < 1000
counter = counter + 1;    
N = 500; % define sequence length
Index = 'AGCT';
rand_seq = '';
count = 0
while count < 500
    Vari = Index(randi(numel(Index)));
    rand_seq = strcat(rand_seq, Vari);
    count = numel(rand_seq);
end
SeqSt = strfind(rand_seq, 'ATG');
SeqEd1 = strfind(rand_seq, 'TAA');
SeqEd2 = strfind(rand_seq, 'TGA');
SeqEd3 = strfind(rand_seq, 'TAG');
SeqEdAll = [SeqEd1, SeqEd2, SeqEd3];
SeqEd = sort(SeqEdAll);
x = numel(SeqSt);
b = numel(SeqEd);
y = 1;
SeqNew = [];
SeqFinal = [];
while y<(x+1)
    a = 1;
    while a<(b+1)
    Sub = SeqSt(y)- SeqEd(a);
    SeqNew = [SeqNew, Sub];
    a = a + 1;
    end
    Seq2 = SeqNew(SeqNew>=0);
    Seq2 = sort(Seq2);
    if numel(Seq2) < 1
    SeqFinal = SeqFinal
    else
    SeqFinal = [SeqFinal, Seq2(1)]
    end
    y = y + 1;
end
SeqFinal = fliplr(sort(SeqFinal));
if numel(SeqFinal) < 1
LongORF = 'undefined'
else
LongORF = SeqFinal(1);
end
%%Probability Array Maker
if LongORF >= 50
    ProbArray = [ProbArray, 1];
else if LongORF == 'undefined'
        ProbArray = [ProbArray, 0];
    else ProbArray = [ProbArray, -1];
    end
end
end
FinalProb = numel(ProbArray(ProbArray == 1))/1000*100
%%My answer is 3.7% for the final probability.
%%
%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 
N = 50;
PlotX = [];
PlotY = [];
while N<500 % Repeat iteration by bp number increments by 50
counter = 1;
ProbArray = [];
while counter < 1000 %Repeat 1000 times to calculate average probability of a >= 50bp ORF
counter = counter + 1; 
Index = 'AGCT';
rand_seq = '';
count = 0;
while count < N %creation of the random sequence of variable length
    Vari = Index(randi(numel(Index)));
    rand_seq = strcat(rand_seq, Vari);
    count = numel(rand_seq);
end
SeqSt = strfind(rand_seq, 'ATG'); %defining arrays used to extrapolate the largest span of ORF
SeqEd1 = strfind(rand_seq, 'TAA');
SeqEd2 = strfind(rand_seq, 'TGA');
SeqEd3 = strfind(rand_seq, 'TAG');
SeqEdAll = [SeqEd1, SeqEd2, SeqEd3];
SeqEd = sort(SeqEdAll);
x = numel(SeqSt);
b = numel(SeqEd);
y = 1;
SeqNew = [];
SeqFinal = [];
while y<(x+1) %calculates each ATG vs the stop codons and takes the smallest positive value to insert into array
    a = 1;
    while a<(b+1)
    Sub = SeqSt(y)- SeqEd(a);
    SeqNew = [SeqNew, Sub];
    a = a + 1;
    end
    Seq2 = SeqNew(SeqNew>=0);
    Seq2 = sort(Seq2);
    if numel(Seq2) < 1
    SeqFinal = SeqFinal;
    else
    SeqFinal = [SeqFinal, Seq2(1)];
    end
    y = y + 1;
end
SeqFinal = fliplr(sort(SeqFinal));
if numel(SeqFinal) < 1
LongORF = 0;
else
LongORF = SeqFinal(1); %after full Start vs Stop array comparison, this longest ORF value is inserted here.
end
%%Probability Array Maker
if LongORF >= 50
    ProbArray = [ProbArray, 1];
else if LongORF == 0
        ProbArray = [ProbArray, 0];
    else ProbArray = [ProbArray, -1];
    end
end
end
FinalProb = numel(ProbArray(ProbArray == 1))/1000*100 %calculates final probability
PlotX = [PlotX, N];
PlotY = [PlotY, FinalProb];
N = N + 50;
end
plot(PlotX,PlotY)%plots the final probability, see values of PlotY for the curve per 50 base pairs.
%appears to peak around 250 bp or so and plateaus

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

%See above for my plot explanation. Appears to be correct as statistically
%after a certain length, the probability of >50 bps is limited by the
%occurence of stop codon vs start codon probility and not by the limited
%length (~250 from my plot).

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

cd 'C:\Users\Superjus91\Documents\GitHub\hw1-jabatin'; %change to directory where qPCRdata.txt exists.
fid = fopen('qPCRdata.txt');
textData = textscan(fid, '%s%s%s%s%s%s%s%s',72, 'headerlines', 2, 'CommentStyle', '> -','EmptyValue',Inf);
%Tab Delimiter does not appear to work, but it does if I replace tab with
%commas in the original file... this causes some wells to be mislabled.
fclose(fid);
CVArray1 = [textData{6}];
CVArray2 = [textData{3},textData{6}]; %place Pos and Cp into an 2x74 Array.

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 

CVArrayR = transpose(reshape(CVArray1,[12,6]));

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

%determine the dimensions of array
CVArrayD = str2double(CVArrayR);
Arow = size(CVArrayR, 1);
Acolumn = size(CVArrayR, 2);
Rep = 3 %specify number of replicates
dCtA = [];
dStd = [];
num2 = 1; %obtain average of the 3 triplicates into a new array
while num2 <= Arow
    num = 0;
while num < Acolumn 
    dCtA = [dCtA, ((CVArrayD(num2, num+1)+CVArrayD(num2, num+2)+CVArrayD(num2, num+3))/Rep)]; 
    dStd = [dStd, std([CVArrayD(num2, num+1),CVArrayD(num2, num+2),CVArrayD(num2, num+3)])];
    num = num + Rep;
end
    num2 = num2 + 1;
end
AvgdCt = transpose(reshape(dCtA, [4, 6])); %the new averaged array
dStdN = transpose(reshape(dStd, [4, 6])); %repeat for Std Dev
AvgdCtN = AvgdCt;
ArowN = 6;
num2 = 1;
while num2 <= ArowN
    num = 1;
while num <= 4 
    AvgdCtN(num2, num) = AvgdCt(num2, num) - AvgdCt(num2, 4);  
    num = num + 1;
end
    num2 = num2 + 1;
end

ddCt = reshape(AvgdCtN(AvgdCtN ~= 0), [6,3]); %now the last row subtract out the control (column 1).
num2 = 1; %obtain average of the 3 triplicates into a new array
while num2 <= ArowN
    num = 1;
while num <= 3 
    ddCt(num2, num) = AvgdCtN(num2, num) - AvgdCtN(num2, 1);  
    num = num + 1;
end
    num2 = num2 + 1;
end
ddCtN = reshape(ddCt(ddCt ~= 0), [6,2]); %remove first row
ddCtNN = ddCtN;
% calculate fold change
num2 = 1;
while num2 <= ArowN
    num = 1;
while num <= 2 
    ddCtNN(num2, num) = 2^-(ddCtN(num2, num));  
    num = num + 1;
end
    num2 = num2 + 1;
end
bar(ddCtNN)
xlabel('Genes')
ylabel('Fold Change')
legend('Treatment1', 'Treatment2', 'Location', 'Northwest')

%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?


% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty

%Answer, added in calcs for Std Dev, but did not incorporate into graph.
%Solutions online suggest overlaying a separate graph with same dimensions.


