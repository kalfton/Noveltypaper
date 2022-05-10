# Novelty paper
This is the repository for sharing the code for the paper: [*"Surprise and recency in novelty detection in the
primate brain"*](https://www.cell.com/current-biology/pdfExtended/S0960-9822(22)00504-8) by Kaining Zhang, Ethan S. Bromberg-Martin, Fatih Sogukpinar, Kim Kocher and Ilya E. Monosov.


To run the analyses and generate the figures in the paper, just run the MATLAB script "Mega_script.m".
The output result plots are saved under the folder "plots".


To plot most of the analysis in the paper, the data in the "maindata" folder is enough. The two .mat files each contain a struct array (Neuronlist_\*) for each animal, and in the struct array, each item contains the precalculated results of one neuron.

Key fields in the Neuronlist: 

| Field name | Description |
| --- | ----------- |
|'name' | Neuron name |
|'recency_ind_match_pos' | Recency index |
|'P_recency_ind_match_pos' | P-value of recency index |
|'pred_nov_vs_fam' | Novelty index |
|'P_pred_nov_vs_fam' | P-value of novelty index |
|'violation_ind' | Violation index |
|'P_violation_ind_perm' | P-value of violation index |
|'learning' | Storing the firing rates to the fractal objects in the learning trial type (z-scored) |
|'nonlearning' | Storing the firing rates to the fractal objects in the non-learning trial types (z-scored) |
|'All_Fractal_FR' | Storing the firing rates to all objects (raw firing rate) |
|'filename' | Name of the raw file that the neuron is in. |
|'zscorestd' | Z-scored standard deviation |
|'zscoremean' | Z-scored mean |
|'monkeyName'  | Monkey name |
|'MonkeyID' | Monkey ID |
|'region'  | The brain area where the neuron belongs. |
|'regionIndex' | Index of the region |
|'electrodeLocation' | Electrode location [lateral, anterior, dorsal] in millimeter |
|'rewardvalueindex_precue' | Reward index |
|'rewardvalueindexP_precue' | P-value of the reward index |
|'RewInfoAnticipIndex_split' | Information anticipation index |
|'RewInfoAnticipIndexP_split' | P-value of the information anticipation index |
|'Learning_\*_start' | Averaged firing rate at the start of the session.(z-scored) |
|'Learning_\*_end' | Averaged firing rate at the end of the session. (z-scored) |
|'Learning_\*_startroc' | Area under ROC of the firing rates to the learning/novel fractals with the firing rate to the familiar fractals at the start of the session. |
|'Learning_\*_endroc' | Area under ROC of the firing rates to the learning/novel fractals with the firing rate to the familiar fractals at the end of the session. |
|'Learning_\*_startp' | P-value of Area under ROC at the start of the session. |
|'Learning_\*_endp' | P-value of Area under ROC at the end of the session. |
|'learning_surprise' | Surprise index in learning, related to figure S6A |
|'withindaylearning_newindex' | Within-day learning index |
|'acrossdayforget_newindex' | Across-day forgetting index |
|'withindaylearning_newindex_p' | P-value of with-inday learning index |
|'acrossdayforget_newindex_p' | P-value of across-day forgetting index |
|'learningforgetinganalysis'  | Criterion to be included in the learning-forgetting correlation analysis |

</br>

## Raw data

The raw neuronal spike data can be downloaded [here](https://wustl.box.com/s/v4x3zjvyopexyud3ghnk87apav9ma6ay).

Each raw file (task_session_\*.mat) is a recording session, and it contains several struct variables, The struct variables whose names look like 'SPK*' contain the sorted neuron's information in each session. The information includes the neuron's spike time point in each trial (field name: sptimes) and each fractal object (field name: fracsptimes). 

The raw file also contains a task information struct named 'Generaltask'. Generaltask contains the behavior procedure information.

Key fields in the Generaltask:

| Field name | Description |
| --- | ----------- |
|'trialstart' | Trial start time |
|'trialend' | Trial end time |
|'successtrial' | Whether the trial had been successfully completed |
|'trialtype' | Trial type. Numbered from 1 to 6: 1 corresponds to trial type 4 in the paper, 2 & 3 correspond to trial type 2, 4 & 5 correspond to trial type 3, and 6 corresponds to trial type 1.
|'trialnumber' | Trial number |
|'timetargeton' | The time that the fractal objects appeared, it is an n by 3 matrix, n is the number of trials,  each row represents a trial and has the 3 fractal objects appearing times. The trial start time has been subtracted from this time|
|'Set' | The fractal object ID in each trial. 7999: novel fractal objects; 63xx, 64xx, 65xx: familiar fractal objects; 73xx, 74xx: repeating novel objects. |
|'timetargetoff' | The time that the fractal objects disappeared, it is an n by 3 matrix, n is the number of trials,  each row represents a trial and has the 3 fractal objects disappearing times. The trial start time has been subtracted from this time |
|'timefpon' | The time that the fixation point appeared, each row represents a trial, the first column represents the fixation point in the object sequence viewing part, and the second column represents the reward fixation point in the instrumental behavior part. The trial start time has been subtracted from this time |
|'timefpoff' | The time that the fixation point appeared. Each row represents a trial, the first column represents the fixation point in the object sequence viewing part, and the second column represents the reward fixation point in the instrumental behavior part. The trial start time has been subtracted from this time |
|'fixacq' | The time that the monkey started to fixate at the fixation point. The trial start time has been subtracted from this time |
|'b_error_place' | Whether the sequence violation happened, and where it happened. 0 represents no sequence violation occurred, 2 represents that the sequence violation occurred in the second place, 3 represents that sequence violation occurred in the third place. |
|'Fractals' | It is a 5 by m matrix, where m is the total number of fractal objects in the session. The first row is fractal ID. The second row is the trial type that the fractal is in, with the same convention as 'trialtype' field. The third row is the fractal position in the sequence. The fourth row is the fractal appearing time. The fifth row is the object disappearing time. |

The raw data will be needed to regenerate some intermediate results, in that case, they should be put in the "Plot_all_session_code/raw_data/" directory.

