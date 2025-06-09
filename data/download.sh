root_url="https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.2"

curl -O "${root_url}/DMS_ProteinGym_substitutions.zip"
unzip -j DMS_ProteinGym_substitutions.zip -d substitutions
rm DMS_ProteinGym_substitutions.zip

curl -O "${root_url}/ProteinGym_AF2_structures.zip"
unzip -j ProteinGym_AF2_structures.zip -d structures
rm ProteinGym_AF2_structures.zip

curl -O "${root_url}/zero_shot_substitutions_scores.zip"
unzip zero_shot_substitutions_scores.zip -d scores
mv scores/zero_shot_substitutions_scores/* scores
rm -r scores/zero_shot_substitutions_scores/
rm zero_shot_substitutions_scores.zip

curl "${root_url}/DMS_substitutions.csv" > dms.csv
