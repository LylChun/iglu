# (1) apply to each subject (episode_single)
#       (a) interpolate to create equidistant grid
#       (b) categorize each data point & label contiguous segments split by gaps
#       (c) classify each segment (event_class)
#       (d) summarise episodes (episode_summary)
# (a) event_class: label events of each type for each segment
#       (a) must be >= duration (function input is # idx to get # minutes)
#       (b) ends at >= dur_length (function input is top level dur_length/dt0)
# (b) episode_summary: calculate summary statistics
#       (a) return for each type of episode: # episodes, mean duration, mean glu value


# helpers
vseq = Vectorize(seq.default, c("from", "to"))

# function to classify and label all events in a segment
event_class = function(data, level_label, event_duration, end_duration) {

  # define ends for each type
  end_cats = switch(
    level_label,
    "lv1_hypo" = c("normal", "lv1_hyper", "lv2_hyper"),
    "lv2_hypo" = c("lv1_hypo", "normal", "lv1_hyper", "lv2_hyper"),
    "lv1_hyper" = c("normal", "lv1_hypo", "lv2_hypo"),
    "lv2_hyper" = c("lv1_hyper","normal", "lv1_hypo", "lv2_hypo")
  )

  # group to include lv2 as a subset of lv1
  # i.e. if lv1 event, summarize to (lv1 U lv2) or not
  # and if lv2 event, summarize to lv2 or not
  level_group = dplyr::if_else(
    rep(grepl("lv1", level_label), nrow(data)),
    dplyr::if_else(data$level %in% end_cats, "end",
                   sub("lv[1|2]_", "", level_label, perl = TRUE)),
    dplyr::if_else(data$level %in% end_cats, "end", data$level)
  )
  level_label = dplyr::if_else(
    grepl("lv1", level_label),
    sub("lv[1|2]_", "", level_label, perl = TRUE),
    level_label
  )

  annotated = data %>%
    # summarize as end type or not
    dplyr::mutate(
      level = level_group
    ) %>%
    # create grouping for each "event" with grouping defined above
    dplyr::mutate(
      event = rep(1:length(rle(level)$lengths), rle(level)$lengths)
    ) %>%
    dplyr::group_by(event) %>%
    dplyr::mutate(
      # possibly category (incl level) event; where duration depends on lv# or extended
      pos_start = any(level == level_label) && (dplyr::n() >= event_duration),
      # if possible event, add start on first index of event
      start = c(dplyr::if_else(pos_start[1], "start", NA_character_),
                rep(NA_character_, dplyr::n()-1)),
      # add possible ends (always need to check for end duration)
      pos_end = any(level == "end") && (dplyr::n() >= end_duration),
      end = c(dplyr::if_else(pos_end[1], "end", NA_character_),
              rep(NA_character_, dplyr::n()-1))
    ) %>%
    dplyr::ungroup()

  ### for each possible end find the matching start
  starts = which(!is.na(annotated$start))
  # add 0 for initial "previous end" to be 0 and length + 1 for last point to be end
  ends = c(0, which(!is.na(annotated$end)), nrow(data) + 1)

  # empty table for loop
  pairs = tibble::tibble(
    start = rep(NA, length(ends)),
    end = rep(NA, length(ends))
  )

  # match end to the minimum start before it that is not linked to a previous end
  for (i in 1:length(ends) - 1) {
    end_idx = i + 1

    # start must occur before this end AND after previous end
    event_bool = starts < ends[end_idx] & starts > ends[end_idx - 1]
    if (any(event_bool)) {
      # choose first start: if multiple starts interrupted but w/o end, choose first start
      pairs$start[i] = starts[min(which(event_bool))]
      # end as last measurement before end
      # (i.e. 60, 65, 75 end on 65 which is idx before normal)
      pairs$end[i] = ends[end_idx] - 1
    }
  }

  # drop nas and create sequence of indices
  pairs = pairs[complete.cases(pairs), ]

  # no episodes in this segment
  if (nrow(pairs) == 0) {
    output = rep(0, nrow(data))
    return (output)

  } else {

    event_idx = unlist(vseq(pairs$start, pairs$end))
    event_label = rep(1:nrow(pairs), (pairs$end - pairs$start)+1)

    output = rep(0, nrow(data))
    output[event_idx] = event_label

    # return event vector with length = nrow(data); not event = 0, event given group label
    # label is unique per subject-segment
    return(output)
  }
}

# function to calculate summary statistics for each type of episode
episode_summary = function (data, dt0) {

  episode_summary_helper = function (data, level_label, dt0) {

    data = data[, c(1:5, which(colnames(data) == level_label))]
    colnames(data) = c("id", "time", "gl", "segment", "level", "event")

    # if none of this event exists, return NA and exit
    if (all(data$event == 0)) {
      output = tibble::tibble(
        avg_ep_per_day = NA_real_,
        avg_ep_duration = NA_real_,
        avg_ep_gl = NA_real_
      )
      return(output)
    }

    data_sum = data %>%
      dplyr::filter(event != 0) %>%
      dplyr::group_by(segment, event) %>%
      dplyr::summarise(
        # for each event, pull duration in minutes and mean glucose
        event_duration = dplyr::n() * dt0,
        event_glucose = mean(gl),
        .groups = "drop"
      )  %>%
      # summarise across all segments and all events
      dplyr::summarise(
        total_events = dplyr::n(),
        total_time = nrow(data) * dt0, # in minutes
        # divide total number of episodes by recording time, not including gap times
        avg_ep_per_day = total_events/(total_time/60/24),
        avg_ep_duration = mean(event_duration),
        avg_ep_gl = mean(event_glucose)
      )

    output = as.vector(data_sum[, 3:5])
    return(output)
  }


  output = tibble::tibble(
    type = c(rep("hypo", 3), rep("hyper", 3)), # hypo or hyper
    level = rep(c("lv1", "lv2", "extended"), 2), # lv1/lv2/extended
    avg_ep_per_day = NA_real_,
    avg_ep_duration = NA_real_,
    avg_ep_gl = NA_real_
  )

  labels = c("lv1_hypo", "lv2_hypo", "ext_hypo", "lv1_hyper", "lv2_hyper", "ext_hyper")
  for (i in 1:length(labels)) {
    output[i, 3:5] = episode_summary_helper(data, labels[i], dt0)
  }

  return(output)
}

# classify episodes for all segments for one subject
episode_single = function(data, lv1_hypo, lv2_hypo, lv1_hyper, lv2_hyper,
                          dur_length, dt0, inter_gap, tz) {

  ### interpolate and segment to deal with gaps and uneven grid
  data_ip <- CGMS2DayByDay(data, dt0 = dt0, inter_gap = inter_gap, tz = tz)
  # find first day and number of days
  day_one = lubridate::as_datetime(data_ip[[2]][1], tz = tz)
  ndays = length(data_ip[[2]])
  # generate grid times by starting from day one and cumulatively summing
  time_ip =  day_one + lubridate::minutes(cumsum(
    # replicate dt0 by number of measurements (total minutes/dt0)
    rep(data_ip$dt0, ndays * 24 * 60 /data_ip$dt0)))

  # interpolated dataframe
  new_data = tibble::tibble(
    id = data$id[1],
    time = time_ip,
    # t to get rowwise vector of a matrix
    gl = as.vector(t(data_ip$gd2d))
  )

  dur_idx = dur_length/data_ip$dt0

  ### add type and level label
  segment_data <- new_data %>%
    dplyr::mutate(
      na = is.na(gl),
      segment = rep(1:length(rle(na)$lengths), rle(na)$lengths)
    ) %>%
    dplyr::filter(!na) %>%
    dplyr::select(-na) %>%
    dplyr::mutate(
      level = dplyr::case_when(
        # below level 2 hypoglycemia
        gl < lv2_hypo ~ "lv2_hypo",
        # in level 1 hypoglycemia (between lv1 and lv2 threshold)
        lv2_hypo <= gl & gl < lv1_hypo ~ "lv1_hypo",
        # above level 2 hyperglycemia
        gl > lv2_hyper ~ "lv2_hyper",
        # in level 1 hyperglycemia (between lv1 an lv2 threshold)
        lv1_hyper < gl & gl <= lv2_hyper ~ "lv1_hyper",
        # normal
        TRUE ~ "normal"
      )
    )

  ep_per_seg = segment_data %>%
    dplyr::group_by(segment) %>%
    dplyr::mutate(
      lv1_hypo = event_class(data.frame(id, time, gl, level), "lv1_hypo", dur_idx, dur_idx),
      lv2_hypo = event_class(data.frame(id, time, gl, level), "lv2_hypo", dur_idx, dur_idx),
      lv1_hyper = event_class(data.frame(id, time, gl, level), "lv1_hyper", dur_idx, dur_idx),
      lv2_hyper = event_class(data.frame(id, time, gl, level), "lv2_hyper", dur_idx, dur_idx),
      # extended hypoglycemia defined as >120 minutes of hypoglycemia
      ext_hypo = event_class(data.frame(id, time, gl, level), "lv1_hypo", 120/data_ip$dt0 + 1, dur_idx),
      ext_hyper = 0 # placeholder for now
    )

  output = episode_summary(ep_per_seg, data_ip$dt0)

  return(output)
}


# top level function overall
episode_refactored = function(data, lv1_hypo = 100,lv2_hypo = 70, lv1_hyper= 120, lv2_hyper = 180,
                              dur_length = 15, dt0 = NULL, inter_gap = 45, tz = "") {

  if (dur_length > inter_gap) {
    warning("Interpolation gap parameter less than episode duration, data gaps may cause
            incorrect computation")
  }

  out <- data %>%
    dplyr::group_by(id) %>%
    # calculate for each subject
    dplyr::reframe(episode_single(data.frame(id, time, gl), lv1_hypo, lv2_hypo, lv1_hyper,
                                  lv2_hyper, dur_length, dt0, inter_gap, tz)) %>%
    dplyr::ungroup()

  return(out)

}
