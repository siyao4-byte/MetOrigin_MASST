# Graph 1 from complete indexed evidence; display grouping never changes exports.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
run_dir <- Sys.getenv("METORIGIN_RUN_DIR", ".")
td <- file.path(run_dir, "outputs/metorigin_matches")
od <- file.path(run_dir, "Graph1_MetOrigin_R_outputs")
dir.create(od, recursive = TRUE, showWarnings = FALSE)
associations <- read_csv(file.path(td, "source_associations.csv"), show_col_types = FALSE)
positive <- associations %>% filter(status == "Yes")
lineages <- read_csv(file.path(td, "bacterial_lineages.csv"), show_col_types = FALSE)
bacterial <- lineages %>% transmute(feature_id, identity_id, compound_name, record_id,
                                  source = "bacteria", detail = organism, layer3 = coalesce(phylum, "Unspecified")) %>% distinct()
other <- positive %>% filter(source != "bacteria") %>%
  mutate(layer3 = coalesce(detail, "Unspecified")) %>%
  select(feature_id, identity_id, compound_name, record_id, source, detail, layer3)
# Malformed taxonomy retains its organism association, with no guessed lineage.
unresolved <- positive %>% filter(source == "bacteria") %>%
  anti_join(bacterial, by = c("feature_id", "record_id", "detail")) %>%
  mutate(layer3 = "Taxonomy unresolved") %>%
  select(feature_id, identity_id, compound_name, record_id, source, detail, layer3)
long <- bind_rows(other, bacterial, unresolved) %>%
  mutate(general_class = str_to_title(source), feature_label = compound_name,
         layer4 = coalesce(detail, "Unspecified"), layer3_plot = layer3, layer4_plot = layer4) %>% distinct()
write_csv(long, file.path(od, "Graph1a_feature_general_subclass_detail_long_table.csv"))
members <- positive %>% distinct(identity_id, source)
sources <- c("human", "bacteria", "food", "plant", "animal", "drug", "environment")
colours <- c(human="#C77979", bacteria="#5B9FC6", food="#D4B85E", plant="#79A774",
             animal="#AB8C73", drug="#9B80BF", environment="#85969B")
theme_set(theme_minimal(base_size = 12) + theme(plot.title = element_text(face="bold"),
                                              panel.grid.minor=element_blank()))
save_plot <- function(p, name, w, h) {
  ggsave(file.path(od, paste0(name,".png")), p, width=w, height=h, dpi=300, bg="white")
  ggsave(file.path(od, paste0(name,".pdf")), p, width=w, height=h, device=cairo_pdf, bg="white")
}
# Every identity and all of its general-source edges remain represented.
nodes <- members %>% group_by(identity_id) %>% summarise(primary=first(source), .groups="drop") %>%
  arrange(match(primary, sources), identity_id) %>% mutate(theta=2*pi*(row_number()-1)/n(), x=cos(theta),y=sin(theta))
centres <- tibble(source=sources,theta=2*pi*(seq_along(sources)-1)/length(sources),
                  x0=.40*cos(theta),y0=.40*sin(theta)) %>%
  left_join(members %>% count(source), by="source") %>% mutate(n=replace_na(n,0L))
edges <- members %>% left_join(nodes,by="identity_id") %>% left_join(centres,by="source")
circular <- ggplot() + geom_segment(data=edges,aes(x=x,y=y,xend=x0,yend=y0,colour=source),alpha=.035,linewidth=.12) +
  geom_point(data=nodes,aes(x,y,colour=primary),size=.25) +
  geom_point(data=centres,aes(x0,y0,fill=source,size=n),shape=21,colour="white") +
  geom_text(data=centres,aes(x0,y0+.10,label=paste0(str_to_title(source),"\n",n)),size=3.5) +
  scale_size_area(max_size=17) + scale_colour_manual(values=colours) + scale_fill_manual(values=colours) +
  coord_equal(xlim=c(-1.12,1.12),ylim=c(-1.12,1.12)) + theme_void() + theme(legend.position="none", plot.title=element_text(face="bold")) +
  labs(title="Recorded source associations")
save_plot(circular,"Graph1a_2layer_circular_center_nodes",10,10)

sets <- members %>% group_by(identity_id) %>% summarise(combination=paste(sources[sources %in% source],collapse=";"),.groups="drop")
counts <- sets %>% count(combination,sort=TRUE) %>% mutate(intersection=row_number())
write_csv(counts,file.path(od,"Graph1b_all_intersections.csv"))
top <- slice_head(counts,n=25)
matrix <- crossing(intersection=top$intersection,source=sources) %>% left_join(top,by="intersection") %>%
  rowwise() %>% mutate(present=source %in% str_split(combination,";",simplify=FALSE)[[1]]) %>% ungroup()
bars <- ggplot(top,aes(intersection,n)) + geom_col(fill="#60666A") +
  geom_text(aes(label=n),vjust=-.4,size=3) + scale_y_continuous(expand=expansion(mult=c(0,.14))) +
  scale_x_continuous(limits=c(.5,nrow(top)+.5),expand=expansion(mult=0)) +
  labs(x=NULL,y="Unique identities") + theme(axis.text.x=element_blank(),panel.grid.major.x=element_blank())
dots <- ggplot(matrix,aes(intersection,factor(source,levels=rev(sources)))) +
  geom_line(data=filter(matrix,present),aes(group=intersection),linewidth=.4) +
  geom_point(colour="#DDDDDD",size=2.5) + geom_point(data=filter(matrix,present),size=2.5) +
  scale_x_continuous(limits=c(.5,nrow(top)+.5),expand=expansion(mult=0)) +
  labs(x="Recorded-source intersections (top 25)",y=NULL) + theme(panel.grid=element_blank(),axis.text.x=element_blank())
save_plot(bars/dots + plot_layout(heights=c(3,1.5)),"Graph1b_origin_upset",14,8)

# This is an association network: links are not a conserved flow of unique molecules.
detail_counts <- long %>% distinct(identity_id,source,layer3) %>% count(source,layer3,name="n")
write_csv(detail_counts,file.path(od,"Graph1a_all_specific_association_counts.csv"))
display <- detail_counts %>% group_by(source) %>% arrange(desc(n),layer3,.by_group=TRUE) %>%
  mutate(label=if_else(row_number()<=8,layer3,"Other recorded subclasses")) %>% ungroup()
grouped <- long %>% distinct(identity_id,source,layer3) %>% left_join(select(display,source,layer3,label),by=c("source","layer3")) %>%
  distinct(identity_id,source,label) %>% count(source,label,name="n") %>%
  arrange(match(source,sources),desc(n)) %>%
  mutate(display_label=str_wrap(paste0(label," (",n,")"),48),
         label_lines=str_count(display_label,"\n")+1L,
         row_height=pmax(1.4,label_lines*.8),
         source_gap=if_else(source != lag(source,default=first(source)),1.2,0),
         y=cumsum(row_height+source_gap))
left <- members %>% count(source,name="total") %>% left_join(grouped %>% group_by(source) %>% summarise(y0=mean(y),.groups="drop"),by="source")
links <- grouped %>% left_join(left,by="source")
identity_nodes <- nodes %>% mutate(y1=seq(min(grouped$y),max(grouped$y),length.out=n()))
incoming <- members %>% left_join(select(identity_nodes,identity_id,y1),by="identity_id") %>% left_join(left,by="source")
sankey <- ggplot(links) +
  geom_segment(data=incoming,aes(x=-1,y=y1,xend=0,yend=y0,colour=source),alpha=.025,linewidth=.12) +
  geom_point(data=identity_nodes,aes(-1,y1),size=.2,colour="#555555") +
  annotate("text",x=-1.12,y=0,label=paste(nrow(identity_nodes),"identities\n(names in source table)"),size=3,hjust=.5) +
  geom_curve(aes(x=0,y=y0,xend=1,yend=y,linewidth=n,colour=source),curvature=.12,alpha=.45) +
  geom_point(aes(1,y,colour=source),size=3) +
  geom_text(aes(1.04,y,label=display_label),hjust=0,size=3) +
  geom_text(data=left,aes(-.04,y0,label=paste0(str_to_title(source),"\n",total)),hjust=1,size=3.5) +
  scale_colour_manual(values=colours) + scale_linewidth_continuous(range=c(.3,7)) +
  scale_x_continuous(limits=c(-1.65,2.6)) + scale_y_reverse(expand=expansion(add=2)) + theme_void() +
  theme(legend.position="none",plot.title=element_text(face="bold")) + labs(title="Recorded source associations")
save_plot(sankey,"Graph1a_3layer_sankey",15,max(9,max(grouped$y)*.23))
writeLines(c("Counts represent unique compound-name identities (SMILES fallback).",
             "Associations overlap: a compound can contribute to multiple sources and body sites.",
             "Only the display groups subclasses beyond the top eight; exports retain all records.",
             "Human-only recorded association does not demonstrate exclusive human production.",
             "Taxonomy is retained only where original parallel fields can be aligned without guessing."),file.path(od,"CAPTION.md"))
