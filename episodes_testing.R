library(iglu)
library(dplyr)

source("R/episode_refactored.R")

new_episodes = episode_refactored(example_data_5_subject)
original_episodes = episode_calculation(example_data_5_subject)

# comparison, first column is original, second column is new
hypo_idx = c(1, 7, 13, 19, 25, 2, 8, 14, 20, 26)
hypo_ep_number = cbind(original_episodes[, 1:2],
                       new_episodes[hypo_idx, 4]) %>%
  mutate(
    difference = avg_ep_per_day - Hypo_ep,
    percent = difference/avg_ep_per_day*100
  )

hyper_idx = c(4, 10, 16, 22, 28, 5, 11, 17, 23, 29)
hyper_ep_number = cbind(original_episodes[, c(1, 3)],
                        new_episodes[hyper_idx, 4]) %>%
  mutate(
    difference = avg_ep_per_day - Hyper_ep,
    percent = difference/avg_ep_per_day*100
  )

hypo_ep_duration = cbind(original_episodes[, c(1, 4)],
                         new_episodes[hypo_idx, 5]) %>%
  mutate(
    difference = avg_ep_duration - hypo_duration,
    percent = difference/avg_ep_duration*100
  )

hyper_ep_duration = cbind(original_episodes[, c(1, 5)],
                          new_episodes[hyper_idx, 5]) %>%
  mutate(
    difference = avg_ep_duration - hyper_duration,
    percent = difference/avg_ep_duration*100
  )



# other function outputs are not comparable, see below

# note we do NOT expect the values to match exactly because the new version:
  # accounts for gaps in the data
  # looks for episodes in contiguous segments not split by calendar day
  # only ends the segment in accordance with the consensus definition
    # i.e. must be >= 15 minutes (dur_length) improvement for episode to end
  # calculates episodes per day counting "days" as collection time, not calendar days
  # additionally adds average glucose by episode type

# further note the new version does not include several of the old outputs:
  # to my understanding cols 6-10 are redundant with the *_percent iglu metrics

# however one drawback is the new function is noticeably slower
