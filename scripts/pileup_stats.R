args = commandArgs(trailingOnly=TRUE)
filename = args[1]
window_size = as.numeric(args[2])

# filename = "/data-weedomics-1/align_fastq/test/reads/test.pileup"
filename = "/data-weedomics-1/align_fastq/test/reads/group_A.pileup"
window_size = as.numeric("1000")

start_time = Sys.time()

file = file(filename, "r")
window = 0
start_position = 0
vec_chr = c()
vec_window = c()
vec_nloci = c()
vec_MEAN_min_depth = c()
vec_MEAN_max_depth = c()
vec_MEAN_mean_depth = c()
vec_MEAN_sdev_depth = c()
vec_SDEV_min_depth = c()
vec_SDEV_max_depth = c()
vec_SDEV_mean_depth = c()
vec_SDEV_sdev_depth = c()

while ( TRUE ) {
    line = readLines(file, n = 1)
    if ( length(line) == 0 ) {
        break
    }
    # print(line)
    line = unlist(strsplit(line, "\t"))
    chr = line[1]
    pos = as.numeric(line[2])
    cov = as.numeric(line[seq(from=4, to=length(line), by=3)])
    cov_min = min(cov, na.rm=TRUE)
    cov_max = max(cov, na.rm=TRUE)
    cov_mean = mean(cov, na.rm=TRUE)
    cov_sdev = sd(cov, na.rm=TRUE)
    cov_sum = sum(cov, na.rm=TRUE)
    if (is.na(cov_sdev)) {
        cov_sdev = 0.0
    }
    if (window == 0) {
        prev_chr = ""
    } else {
        prev_chr = chr
    }
    if ((pos < start_position+window_size) & (chr == prev_chr)) {
        tmp_min_depth = c(tmp_min_depth, cov_min)
        tmp_max_depth = c(tmp_max_depth, cov_max)
        tmp_mean_depth = c(tmp_mean_depth, cov_mean)
        tmp_sdev_depth = c(tmp_sdev_depth, cov_sdev)
        vec_nloci[window] = vec_nloci[window] + 1
    } else {
        window = window + 1
        vec_chr = c(vec_chr, chr)
        vec_window = c(vec_window, window)
        vec_nloci = c(vec_nloci, 1)
        if (window==1) {
            tmp_min_depth = c(cov_min)
            tmp_max_depth = c(cov_max)
            tmp_mean_depth = c(cov_mean)
            tmp_sdev_depth = c(cov_sdev)
        } else {
            vec_MEAN_min_depth = c(vec_MEAN_min_depth, mean(tmp_min_depth))
            vec_MEAN_max_depth = c(vec_MEAN_max_depth, mean(tmp_max_depth))
            vec_MEAN_mean_depth = c(vec_MEAN_mean_depth, mean(tmp_mean_depth))
            vec_MEAN_sdev_depth = c(vec_MEAN_sdev_depth, mean(tmp_sdev_depth))
            vec_SDEV_min_depth = c(vec_SDEV_min_depth, sd(tmp_min_depth))
            vec_SDEV_max_depth = c(vec_SDEV_max_depth, sd(tmp_max_depth))
            vec_SDEV_mean_depth = c(vec_SDEV_mean_depth, sd(tmp_mean_depth))
            vec_SDEV_sdev_depth = c(vec_SDEV_sdev_depth, sd(tmp_sdev_depth))
        }
        start_position = pos
    }
}
close(file)

end_time = Sys.time()

end_time - start_time