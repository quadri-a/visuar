# visuar
Vital Signs-based User Authentication using mmWave Radar (ViSUAR)

1) Please start by runnning the script user_check.m

2) user_check.m will ask you to enter a case number.
   The database is already populated with 3 cases.
   These cases contain unprocessed heartbeat data collected from 3 users (1 female, 2 male particiapnts).

ViSUAR Flowchart is show below.
![visuar_flowchart](https://github.com/user-attachments/assets/e8cea274-f619-4400-99cf-9187724d8b8d)

Vital signs detection alogorithm include below processes:
1) User's chest motion detection.
2) Phase estimation of received signal provides information on user's heartbeat/breathing related motion.
3) Empirical mode decomposition-based first stage filtering of signal.
4) Second stage filtering includes IIR filtering, which leaves the filtered signal with only heartbeat related data.
5) Cardiac cycle identification and extraction process begins after filtering signal.
6) Feature extraction from cardiac cycle is performed.
7) Principal Component Analysis (PCA) transformation applied on extracted features to find unique heartbeat pattern.
8) Unique heartbeat pattern is stored as user's cardiac profile.
9) Learning algorithm is trained with collected cardiac profile of all users to perform user authentication.
