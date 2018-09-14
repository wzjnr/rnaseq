#{{{
source("me.fun.r")
studies
#fi = '~/data/genome/B73/v32/t5.gtb'
#gids = read_tsv(fi) %>% distinct(par) %>% pull(par)
fi = file.path('~/data/genome/B73', "v37/t2.tsv")
t_gs = read_tsv(fi, col_types = 'ccccciic') %>% 
    filter(etype == 'exon') %>% 
    group_by(gid, tid) %>% 
    summarise(size = sum(end - beg + 1)) %>%
    group_by(gid) %>%
    summarise(size = max(size))
#}}}

#{{{ li2013
study = 'li2013'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts_pa = c("B73", "Mo17")
gts = sort(unique(th$Genotype))
gts = c(gts_pa, gts[! gts %in% gts_pa])

#{{{ number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti %>% mutate(ndiff = input_reads - surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, surviving, dropped) %>%
    gather(type, nseq, -SampleID)
types = c("surviving", "dropped")
tp = th %>% inner_join(ti2, by = 'SampleID') %>%
    mutate(type = factor(type, levels = types),
           Genotype = factor(Genotype, levels = rev(gts))) %>%
    group_by(Genotype, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup()

p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
#ggsave(p1, filename = fp, width = 10, height = 12)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 10)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(SampleID = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 10)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 10)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
#t_rc = tcw %>% gather(sid, RawReadCount, -gid)
t_rc = tcw

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
cols1 = c('gray80','black','red','seagreen3',pal_d3()(5))
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(0,.15)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = SampleID), size = 2, nudge_x = .001, hjust = 0) +
    geom_text(aes(label = Genotype), size = 2, nudge_x= .015, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = .022, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 15)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(gt = ifelse(Genotype %in% c("B73", "Mo17"), Genotype, "RIL"),
           repli = sprintf("Rep%d", Replicate))
p1 = ggplot(tp, aes(x = PC1, y = PC2, color = gt, label = gt, shape = repli)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ hirsch2014
study = 'hirsch2014'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = unique(th$Genotype)

#{{{ number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi) %>%
    replace_na(list('input_reads'=0, 'input_read_pairs'=0,
                    'forward_only_surviving'=0, 'reverse_only_surviving'=0)) %>%
    mutate(input_read_pairs = ifelse(input_read_pairs == 0, input_reads, input_read_pairs))
ti %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
#
ti2 = ti %>% separate(Sample, c("SampleID", "suf"), sep = "_") %>% select(-suf) %>%
    select(SampleID, surviving, forward_only_surviving, reverse_only_surviving, dropped) %>%
    gather(type, nseq, -SampleID)
types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti2, by = 'SampleID') %>%
    mutate(type = factor(type, levels = types),
           Genotype = factor(Genotype, levels = rev(gts))) %>%
    group_by(Genotype, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup()

p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
#ggsave(p1, filename = fp, width = 10, height = 12)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 30)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped) 
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)

types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 30)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 30)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
#t_rc = tcw %>% gather(sid, RawReadCount, -gid)
t_rc = tcw

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
cols1 = c('gray80','black','red','seagreen3',pal_d3()(5))
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits = c(0,1.8)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = Genotype), size = 2, nudge_x= .02, hjust = 0) 
    #geom_text(aes(label = Rep), size = 2, nudge_x = .03, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 30)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(gt = ifelse(Genotype %in% c("B73", "Mo17"), Genotype, "other"))
p1 = ggplot(tp, aes(x = PC1, y = PC2, color = gt, label = gt)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ leiboff2015
study = 'leiboff2015'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
#fh = file.path(dirw, '01.reads.tsv')
fh = file.path(dirw, '02.reads.corrected.tsv')
th = read_tsv(fh)
gts = unique(th$Genotype)

#{{{ trim - number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    select(SampleID, input_read_pairs, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) 
sum(ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>% pull(ndiff))
ti3 = ti2 %>% select(-input_read_pairs) %>%
    gather(type, nseq, -SampleID)

types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .7) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 20)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 20)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 20)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits = c(-.1,5)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = Genotype), size = 2, nudge_x= .02, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = .5, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 30)
#}}}

#{{{ identify mis-labelled replicate
cls = cutree(hc, h = .01)
tcl = tibble(SampleID = names(cls), grp = as.integer(cls))
th2 = th %>% inner_join(tcl, by = 'SampleID')
th3 = th2 %>% group_by(Genotype) %>% 
    summarise(nrep = length(Replicate), ngrp = length(unique(grp))) %>%
    ungroup() %>%
    filter(nrep > 1, ngrp > 1)
th2 %>% filter(Genotype %in% th3$Genotype) %>% print(n=40)
# create 02.reads.corrected.tsv
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)

tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Genotype)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    #scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ jin2016
study = 'jin2016'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = unique(th$Genotype)

#{{{ trim - number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    select(SampleID, input_read_pairs, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) 
sum(ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>% pull(ndiff))
ti3 = ti2 %>% select(-input_read_pairs) %>%
    gather(type, nseq, -SampleID)

types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .7) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 20)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 20)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 20)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits = c(-.2,7.5)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = Genotype), size = 2, nudge_x= .02, hjust = 0) 
    #geom_text(aes(label = Rep), size = 2, nudge_x = .03, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 30)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(gt = ifelse(Genotype %in% c("B73", "Mo17"), Genotype, "other"))
p1 = ggplot(tp, aes(x = PC1, y = PC2, color = gt, label = gt)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ stelpflug2016
study = 'stelpflug2016'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
tissues = sort(unique(th$Tissue))
reps = sort(unique(th$Replicate))

#{{{ number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti %>% mutate(ndiff = input_reads - surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, surviving, dropped) %>%
    gather(type, nseq, -SampleID)
types = c("surviving", "dropped")
tp = th %>% inner_join(ti2, by = 'SampleID') %>%
    mutate(type = factor(type, levels = types),
           Replicate = factor(Replicate, levels = rev(reps)),
           Tissue = factor(Tissue, levels = tissues)) %>%
    group_by(Tissue, Replicate, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup()

p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Replicate, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(SampleID = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Replicate, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Replicate = factor(Replicate, levels = rev(reps)),
           Tissue = factor(Tissue, levels = tissues))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Replicate, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>% ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Replicate, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Replicate = factor(Replicate, levels = rev(reps)),
           Tissue = factor(Tissue, levels = tissues))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Replicate, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>% ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.3,7.5)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = Tissue), size = 3, nudge_x = .02, hjust = 0) +
    geom_text(aes(label = Rep), size = 3, nudge_x = 1.5, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 8)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tismap = LETTERS[1:length(tissues)]
names(tismap) = tissues
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>% 
    mutate(label = tismap[Tissue])
cols = c(brewer.pal(3, 'Set1'))
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Replicate), size = 3) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(name = "", values = tismap) +
    scale_color_manual(name = "", values = cols) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 9, height = 8)
#}}}
#}}}

#{{{ walley2016
study = 'walley2016'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
tissues = sort(unique(th$Tissue))
th = th %>% mutate(Replicate = sprintf("rep%d", Replicate))
reps = sort(unique(th$Replicate))

#{{{ number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti %>% mutate(ndiff = input_reads - surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, surviving, dropped) %>%
    gather(type, nseq, -SampleID)
types = c("surviving", "dropped")
tp = th %>% inner_join(ti2, by = 'SampleID') %>%
    mutate(type = factor(type, levels = types),
           Replicate = factor(Replicate, levels = rev(reps)),
           Tissue = factor(Tissue, levels = tissues)) %>%
    group_by(Tissue, Replicate, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup()

p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Replicate, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(SampleID = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Replicate, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Replicate = factor(Replicate, levels = rev(reps)),
           Tissue = factor(Tissue, levels = tissues))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Replicate, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>% ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Replicate, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Replicate = factor(Replicate, levels = rev(reps)),
           Tissue = factor(Tissue, levels = tissues))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Replicate, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>% ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
#t_rc = tcw %>% gather(sid, RawReadCount, -gid)
t_rc = tcw

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(0,3.1)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = Tissue), size = 3, nudge_x = .02, hjust = 0) +
    geom_text(aes(label = Rep), size = 3, nudge_x = .5, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 8)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tismap = LETTERS[1:length(tissues)]
names(tismap) = tissues
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>% 
    mutate(label = tismap[Tissue])
cols = c(brewer.pal(3, 'Set1'))
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Replicate), size = 3) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(name = "", values = tismap) +
    scale_color_manual(name = "", values = cols) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 9, height = 8)
#}}}
#}}}

#{{{ lin2017
study = 'lin2017'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
#fh = file.path(dirw, '01.reads.tsv')
fh = file.path(dirw, '02.reads.corrected.tsv')
th = read_tsv(fh)
gts = unique(th$Genotype)
tissues = unique(th$Tissue)

#{{{ trim - number read pairs - w. both SE and PE
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% separate(Sample, c("SampleID", "suf"), sep = "_") %>% select(-suf) %>%
    replace_na(list('input_reads'=0, 'input_read_pairs'=0,
                    'forward_only_surviving'=0, 'reverse_only_surviving'=0)) %>%
    mutate(input_read_pairs = ifelse(input_read_pairs == 0, input_reads, input_read_pairs))
ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti3 = ti2 %>% 
    select(SampleID, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) %>%
    gather(type, nseq, -SampleID)

types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .7) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)
#
types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    scale_x_continuous(expand = c(0,0), limits=c(0,15)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) + 
    geom_text(aes(label = SampleID), size = 2, nudge_x = .02, hjust = 0) +
    geom_text(aes(label = Tissue), size = 2, nudge_x = .9, hjust = 0) +
    geom_text(aes(label = Genotype), size = 2, nudge_x = 1.4, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = 2, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 10, height = 15)
# create 02.reads.corrected.tsv
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
gtmap = c(LETTERS[1:26], '@')
names(gtmap) = gts
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>% 
    mutate(label = gtmap[Genotype])
cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
cols = c(brewer.pal(5, 'Set1'))
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, color = Tissue, shape = Genotype), size = 3) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(name = "", values = gtmap) +
    scale_color_manual(name = "", values = cols) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 11, height = 10)
#}}}
#}}}

#{{{ kremling2018
study = 'kremling2018'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = unique(th$Genotype)
tissues = unique(th$Tissue)

#{{{ number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti %>% mutate(ndiff = input_reads - surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, surviving, dropped) %>%
    gather(type, nseq, -SampleID)
types = c("surviving", "dropped")
tp = th %>% inner_join(ti2, by = 'SampleID') %>%
    mutate(type = factor(type, levels = types),
           Genotype = factor(Genotype, levels = rev(gts))) %>%
    group_by(Genotype, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup()

p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
#ggsave(p1, filename = fp, width = 10, height = 12)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 10)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(SampleID = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))

p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 20, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 10)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)
#
types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 20, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=6))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 10)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
#t_rc = tcw %>% gather(sid, RawReadCount, -gid)
t_rc = tcw

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    scale_x_continuous(expand = c(0,0), limits=c(-1,230)) +
    scale_y_discrete(expand = c(.001,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) + 
    geom_text(aes(label = Tissue), size = 2, nudge_x = 1, hjust = 0) +
    geom_text(aes(label = Genotype), size = 2, nudge_x = 11, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = 8, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 10, height = 130, limitsize = F)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') 
cols = c(brewer.pal(7, 'Set1'))

p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, color = Tissue, shape = Tissue), size = 2) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(values = c(15:19,7,8)) +
    scale_color_manual(name = "", values = cols) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 11, height = 10)
#}}}
#}}}

#{{{ baldauf2018
study = 'baldauf2018'

#}}}

#{{{ briggs
study = 'briggs'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = c("B73", "Mo17", "B73xMo17")
tissues = sort(unique(th$Tissue))
th = th %>% filter(Genotype %in% gts, ! SampleID %in% c('BR207', 'BR230', "BR235"))

#{{{ trim - number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    select(SampleID, input_read_pairs, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) 
sum(ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>% pull(ndiff))
ti3 = ti2 %>% select(-input_read_pairs) %>%
    gather(type, nseq, -SampleID)

types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .7) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Genotype = factor(Genotype, levels = rev(gts)),
           Tissue = factor(Tissue, levels = tissues))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 4, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
t_rc = tcw[,c('gid', th$SampleID)]

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)
#
thp = mutate(th, 
    lab = sprintf("%s %s %s %s", SampleID, Tissue, Genotype, Replicate))
p1 = ggtree(tree) + 
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    ggplot2::xlim(0, 16) + 
    theme_tree2()
p1 = p1 %<+% thp + 
    geom_tiplab(aes(label = lab), size = 3, offset = 0.04) + 
    #geom_text(aes(color = as.character(gt_ok), label = gt), size = 4, nudge_x = 6, hjust = 0) + 
    scale_color_manual(values = c("black", "royalblue", "tomato"))
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 12, height = 30)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tismap = LETTERS[1:length(tissues)]
names(tismap) = tissues
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>% 
    mutate(label = tismap[Tissue])
cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
cols = c(brewer.pal(4, 'Set1'))
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Genotype), size = 3) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(name = "", values = tismap) +
    scale_color_manual(name = "", values = cols) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 11, height = 10)
#}}}
#}}}

#{{{ biomAP

#}}}

#{{{ enders 3' RNA-Seq
study = 'enders'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = sort(unique(th$Genotype))
tissues = sort(unique(th$Tissue))

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(SampleID = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Tissue, ncol = 1, scale = 'free') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)
#
types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 1, scale = 'free') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", Tissue, Genotype, Treatment)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.02,3.3)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = lab), size = 2.5, nudge_x = .01, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 10)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
p1 = ggplot(tp, aes(x = PC1, y = PC2, shape = Genotype, color = Tissue)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(1,1), legend.justification = c(1,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ settles endosperm
study = 'settles'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
th = th %>% mutate(Genotype = Treatment)
tissues = sort(unique(th$Tissue))
gts = sort(unique(th$Genotype))

#{{{ trim - number read pairs - w. both SE and PE
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% separate(Sample, c("SampleID", "suf"), sep = "_") %>% select(-suf) %>%
    replace_na(list('input_reads'=0, 'input_read_pairs'=0,
                    'forward_only_surviving'=0, 'reverse_only_surviving'=0)) %>%
    mutate(input_read_pairs = ifelse(input_read_pairs == 0, input_reads, input_read_pairs))
ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti3 = ti2 %>% 
    select(SampleID, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) %>%
    gather(type, nseq, -SampleID)

types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .7) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 1, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Genotype, ncol = 1, scale = 'free') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)
#
types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 1, scale = 'free') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", Tissue, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.02,5.5)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = lab), size = 2.5, nudge_x = .01, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 8)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
p1 = ggplot(tp, aes(x = PC1, y = PC2, shape = Genotype, color = Tissue)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}
#}}}

# combined datasets
#{{{ dev41: stelpflug2016 + walley2016
study = 'dev41'
dirw = file.path(dirp, 'data', study)

#{{{ collect featurecounts data & normalize
studys = c('stelpflug2016', 'walley2016')
th = tibble()
t_rc = tibble()
for (study1 in studys) {
    fh = file.path(dirp, study1, 'data/01.reads.tsv')
    th1 = read_tsv(fh) %>% 
        mutate(Tissue = sprintf("%s_%s", study1, Tissue))
    th = rbind(th, th1)
    fi = file.path(dirp, study1, 'data/raw/featurecounts.tsv')
    t_rc1 = read_tsv(fi)
    if(nrow(t_rc) == 0)
        t_rc = t_rc1
    else {
        if(! identical(t_rc$gid, t_rc1$gid))
            stop('not identical cols[gid]')
        t_rc = t_rc %>% bind_cols(t_rc1[,-1])
    }
}
tissues = sort(unique(th$Tissue))

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fh = file.path(dirw, '01.reads.tsv')
write_tsv(th, fh)
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    scale_x_continuous(expand = c(0,0), limits=c(-.2,12)) +
    scale_y_discrete(expand = c(.02,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) + 
    geom_text(aes(label = Tissue), size = 2, nudge_x = .1, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = 4, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 12)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
cols6 = pal_d3()(6)
shapes7 = c(15:18, 7,8,9)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
tps = tp %>% distinct(Tissue) %>%
    mutate(colors = rep(cols6, each = 7)[1:length(tissues)],
           shapes = rep(shapes7, 6)[1:length(tissues)])
p1 = ggplot(tp) +
    #geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Tissue, group = Tissue), size = 2) +
    geom_point(aes(x = PC1, y = PC2), color = 'red', size = 2) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Tissue), size = 2) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    #scale_shape_manual(breaks = tps$Tissue, values = tps$shapes) +
    #scale_color_manual(breaks = tps$Tissue, values = tps$colors) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 12)
#}}}
#}}}

#{{{ dev64: stelpflug2016 + walley2016 + briggs B73
study = 'dev64'
dirw = file.path(dirp, 'data', study)

#{{{ collect featurecounts data & normalize
studys = c('stelpflug2016', 'walley2016')
th = tibble()
t_rc = tibble()
for (study1 in studys) {
    fh = file.path(dirp, study1, 'data/01.reads.tsv')
    th1 = read_tsv(fh) %>% 
        mutate(Tissue = sprintf("%s_%s", study1, Tissue))
    th = rbind(th, th1)
    fi = file.path(dirp, study1, 'data/raw/featurecounts.tsv')
    t_rc1 = read_tsv(fi)
    if(nrow(t_rc) == 0)
        t_rc = t_rc1
    else {
        if(! identical(t_rc$gid, t_rc1$gid))
            stop('not identical cols[gid]')
        t_rc = t_rc %>% bind_cols(t_rc1[,-1])
    }
}

study1 = 'briggs'
fh = file.path(dirp, study1, 'data/01.reads.tsv')
th1 = read_tsv(fh) %>% 
    filter(Genotype == 'B73', ! SampleID %in% c('BR207', 'BR230', "BR235")) %>%
    mutate(Tissue = sprintf("%s_%s", study1, Tissue))
fi = file.path(dirp, study1, 'data/raw/featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
t_rc1 = tcw[,c('gid', th1$SampleID)]
#
if(! identical(t_rc$gid, t_rc1$gid))
    stop('not identical cols[gid]')
th = rbind(th, th1)
t_rc = t_rc %>% bind_cols(t_rc1[,-1])
tissues = sort(unique(th$Tissue))

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fh = file.path(dirw, '01.reads.tsv')
write_tsv(th, fh)
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    scale_x_continuous(expand = c(0,0), limits=c(-.5,17)) +
    scale_y_discrete(expand = c(.02,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) + 
    geom_text(aes(label = Tissue), size = 2, nudge_x = .1, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = 5, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 18)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
cols6 = pal_d3()(6)
shapes7 = c(15:18, 7,8,9)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
tps = tp %>% distinct(Tissue) %>%
    mutate(colors = rep(cols6, each = 7)[1:length(tissues)],
           shapes = rep(shapes7, 6)[1:length(tissues)])
p1 = ggplot(tp) +
    #geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Tissue, group = Tissue), size = 2) +
    geom_point(aes(x = PC1, y = PC2), color = 'red', size = 2) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Tissue), size = 2) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    #scale_shape_manual(breaks = tps$Tissue, values = tps$shapes) +
    #scale_color_manual(breaks = tps$Tissue, values = tps$colors) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 12)
#}}}
#}}}

