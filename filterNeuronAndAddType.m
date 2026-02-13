% This function splits the sorted neurons into 3 classes based on their
% typical response. This list was created by Pooya by examining each
% neuron. This function is currently required by other sciprts. Tim, 251104

function C = filterNeuronAndAddType(C)
% Build C.neuronType from built-in neuron lists for C.cond.
% Types: 1=Choice–Reward Sustained, 2=Choice-Related Transient, 3=Choice-Suppressed

nSess = numel(C.rawSpikes);
C.neuronType = cell(nSess,1);
for s = 1:nSess
    C.neuronType{s} = zeros(numel(C.rawSpikes{s}),1,'uint8');
end

switch upper(string(C.cond))
    case "LIVE"
        % Choice–Reward Sustained (76)
        sessions1 = [3 3 3 3 3 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
        neurons1  = [2 3 8 12 16 7 8 9 11 4 5 8 9 12 14 2 4 5 6 10 17 3 4 6 8 9 13 15 18 20 21 22 23 29 1 2 3 4 5 7 11 15 16 17 19 22 23 24 27 30 31 32 34 38 39 40 46 47 1 2 3 5 6 10 13 14 15 17 18 21 24 26 29 31 33 36];
        % Choice-Related Transient (34)
        sessions2 = [1 1 1 3 3 3 3 6 6 6 6 7 7 7 7 8 8 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10];
        neurons2  = [1 2 3 1 4 5 24 10 16 17 18 12 13 14 15 7 19 8 9 18 20 21 25 28 43 8 11 20 22 23 25 30 38 39];
        % Choice-Suppressed (9)
        sessions3 = [2 2 6 7 7 7 8 9 9];
        neurons3  = [3 4 3 3 8 9 28 36 45];

    case "REPLAY"
        % Choice–Reward Sustained (59)
        sessions1 = [3 4 4 4 4 4 4 4 4 4 4 4 4 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
        neurons1  = [2 1 2 3 4 5 8 9 10 11 13 15 16 2 3 5 7 15 16 3 4 5 7 9 10 12 15 16 17 18 19 20 21 22 2 3 4 5 6 8 10 11 13 14 15 19 20 21 22 24 25 28 29 30 31 32 33 41 42];
        % Choice-Related Transient (14)
        sessions2 = [5 5 5 5 8 8 8 10 10 10 10 10 10 10];
        neurons2  = [4 8 9 8 1 11 12 7 9 16 27 35 36 38];
        % Choice-Suppressed (12)
        sessions3 = [1 1 5 5 5 5 5 5 5 8 8 10];
        neurons3  = [1 2 1 3 5 6 7 10 11 4 10 40];

    case "AI"
        % Choice–Reward Sustained (70)
        sessions1 = [2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 5 5 5 5 5 5 7 7 7 7 7 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
        neurons1  = [1 3 4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 5 6 9 12 16 17 3 10 15 17 18 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 25 3 4 5 8 9 10 11 12 13 14 16 17 18 19 20 23 24 27 28 33 34];
        % Choice-Related Transient (14)
        sessions2 = [2 2 2 5 5 5 5 7 7 9 10 10 10 10];
        neurons2  = [6 7 8 4 8 11 13 12 14 27 6 7 22 31];
        % Choice-Suppressed (8)
        sessions3 = [2 2 7 7 7 7 9 10];
        neurons3  = [5 9 5 6 8 13 24 32];

    case "DECOY"
        % Choice–Reward Sustained (95)
        sessions1 = [1 1 1 1 2 2 2 2 5 5 5 6 6 6 6 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14];
        neurons1  = [1 2 4 6 1 2 3 5 2 3 5 1 2 3 4 2 5 1 2 3 4 5 7 8 10 5 6 14 15 17 3 5 14 16 17 18 1 2 3 5 6 7 8 9 10 11 12 13 14 15 17 19 2 3 4 5 6 13 15 16 17 19 20 22 24 25 30 31 32 34 35 36 37 39 43 44 3 6 7 10 11 12 13 14 15 19 20 21 22 25 26 27 28 29 31];
        % Choice-Related Transient (27)
        sessions2 = [2 4 8 9 9 9 9 11 11 11 11 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14];
        neurons2  = [9 1 6 13 16 19 20 6 8 12 13 8 9 10 12 18 21 26 27 40 8 9 16 17 18 23 30];
        % Choice-Suppressed (11)
        sessions3 = [1 1 1 2 2 2 4 9 11 13 13];
        neurons3  = [7 9 11 6 8 10 2 11 10 33 42];

    otherwise
        sessions1 = []; neurons1 = [];
        sessions2 = []; neurons2 = [];
        sessions3 = []; neurons3 = [];
end

sessLists = {sessions1, sessions2, sessions3};
neurLists = {neurons1,  neurons2,  neurons3};
for t = 1:3
    ss = sessLists{t}; nn = neurLists{t};
    for i = 1:numel(ss)
        s = ss(i); n = nn(i);
        if s>=1 && s<=nSess && n>=1 && n<=numel(C.rawSpikes{s})
            C.neuronType{s}(n) = uint8(t);
        end
    end
end


end