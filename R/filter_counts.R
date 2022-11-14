# remember to add dplyr to namespace
#' Filter counts based on number of zeros and mean counts, per gene.
#'
#' @param counts Data frame of counts with column `Geneid`
#' @param group Vector of group labels for each sample
#' @param options List with some of `n0`, `n0_grp`, `mean_count`,
#' `mean_count_grp`;
#'
#' @export
filter_counts <- function(counts, group, options) {
    long <- counts %>%
        tidyr::pivot_longer(cols = -c(Geneid),
                            names_to = "sample",
                            values_to = "count")

    # add group column
    group <- rep(group,
                    nrow(counts))
    long$group = group

    # add n0 and mean count columns
    long <- long %>% group_by(Geneid) %>%
        mutate(n0 = sum(count == 0), mean_count = mean(count),
               min_count = min(count)) %>%
        group_by(Geneid, group) %>%
        mutate(n0_grp = sum(count == 0), mean_count_grp = mean(count))

    if (!is.null(options$n0)) {
        long <- long %>% filter(n0 < options$n0)
    }
    if (!is.null(options$mean_count)) {
        long <- long %>% filter(mean_count > options$mean_count)
    }
    if (!is.null(options$n0_grp)) {
        long <- long %>% filter(n0_grp < options$n0_grp)
    }
    if (!is.null(options$mean_count_grp)) {
        long <- long %>% filter(mean_count_grp > options$mean_count_grp)
    }
    if (!is.null(options$min_count)) {
      long <- long %>% filter(min_count >= options$min_count)
    }

    # returns to wide format and
    # drops genes that were only partially filtered by
    # one of the grp filtering options
    df <- long %>% ungroup() %>%
        select(-n0, -mean_count, -group,
               -mean_count_grp, -n0_grp, -min_count) %>%
        pivot_wider(names_from = sample, values_from = count) %>% drop_na()

    df
}
