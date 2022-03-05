# visuar
Vital Signs-based User Authentication using mmWave Radar (ViSUAR)

Please start by runnning the script user_check.m

user_check.m will ask you to enter a case number. The database is already populated with 3 cases. 
These cases contain unprocessed heartbeat data collected from 3 users (1 female, 2 male particiapnts).

Through several steps in the algorithm, chosen user's heartbeat data collected by FMCW radar will be filtered and unique features will be extracted.
For simplicity, in each step a pre-filtered user data for the chosen case is initialized for demonstration.

Vital signs detection alogorithm include below processes:
1) User's chest motion detection (range bi selection).
2) Phase estimation of received signal provides signal that information on user's heartbeat/breathing related motion.
3) Empirical mode decomposition-based first stage filtering of signal.
4) Second stage filtering includes IIR filtering process, which leaves the filtered signal with onyl heartbeat related data.
5) Cardiac cycle identification and extraction process begins after filtering signals.
6) Feature extraction from cardiac cycle starts here.
7) Principal Component Analysis (PCA) transformation applied on extracted features to find unique heartbeat pattern.
8) Unique heartbeat pattern is stored as user's cardiac profile.
9) SVM learning algorithm is trained with collected cardiac profile of all user for user authentication.

Several plots are genreated to demonstrate the result of these processes.
