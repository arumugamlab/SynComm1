#Copy number correction, painful but necessary

#For those with species level resolution, we do the correction manually from the rrndb web page

taxwcpnrs = data.frame(plotpseq_glom@tax_table)

View(taxwcpnrs)

taxwcpnrs$cpnr = NA

taxwcpnrs['Seq27','cpnr'] = 6
taxwcpnrs['Seq3','cpnr'] = 3
taxwcpnrs['Seq383','cpnr'] = 3
taxwcpnrs['Seq165','cpnr'] = 5
taxwcpnrs['Seq88','cpnr'] = 5
taxwcpnrs['Seq67','cpnr'] = 3
taxwcpnrs['Seq64','cpnr'] = 2
taxwcpnrs['Seq86','cpnr'] = 2
taxwcpnrs['Seq247','cpnr'] = 2
taxwcpnrs['Seq353','cpnr'] = 9
taxwcpnrs['Seq259','cpnr'] = 8
taxwcpnrs['Seq2','cpnr'] = 6
taxwcpnrs['Seq5','cpnr'] = 4
taxwcpnrs['Seq7','cpnr'] = 5
taxwcpnrs['Seq8','cpnr'] = 5
taxwcpnrs['Seq9','cpnr'] = 5
taxwcpnrs['Seq22','cpnr'] = 7
taxwcpnrs['Seq29','cpnr'] = 7
taxwcpnrs['Seq31','cpnr'] = 4
taxwcpnrs['Seq37','cpnr'] = 5
taxwcpnrs['Seq39','cpnr'] = 4
taxwcpnrs['Seq59','cpnr'] = 4
taxwcpnrs['Seq75','cpnr'] = 5
taxwcpnrs['Seq87','cpnr'] = 5

taxwcpnrs['Seq28','cpnr'] = 5
taxwcpnrs['Seq80','cpnr'] = 4
taxwcpnrs['Seq152','cpnr'] = 4
taxwcpnrs['Seq191','cpnr'] = 3
taxwcpnrs['Seq201','cpnr'] = 4

taxwcpnrs['Seq180','cpnr'] = 5
taxwcpnrs['Seq65','cpnr'] = 5

taxwcpnrs['Seq77','cpnr'] = 5
taxwcpnrs['Seq11','cpnr'] = 7
taxwcpnrs['Seq73','cpnr'] = 5
taxwcpnrs['Seq173','cpnr'] = 5
taxwcpnrs['Seq97','cpnr'] = 8
taxwcpnrs['Seq114','cpnr'] = 4

taxwcpnrs['Seq26','cpnr'] = 8
taxwcpnrs['Seq53','cpnr'] = 4
taxwcpnrs['Seq76','cpnr'] = 6

taxwcpnrs['Seq171','cpnr'] = 5
taxwcpnrs['Seq107','cpnr'] = 5

taxwcpnrs['Seq395','cpnr'] = 4
taxwcpnrs['Seq150','cpnr'] = 6

taxwcpnrs['Seq1','cpnr'] = 6
taxwcpnrs['Seq122','cpnr'] = 4

taxwcpnrs['Seq23','cpnr'] = 17
taxwcpnrs['Seq24','cpnr'] = 6
taxwcpnrs['Seq49','cpnr'] = 5

taxwcpnrs['Seq403','cpnr'] = 4
taxwcpnrs['Seq352','cpnr'] = 6

taxwcpnrs['Seq4','cpnr'] = 3
taxwcpnrs['Seq125','cpnr'] = 4
taxwcpnrs['Seq14','cpnr'] = 7

taxwcpnrs['Seq92','cpnr'] = 6
taxwcpnrs['Seq108','cpnr'] = 4
taxwcpnrs['Seq21','cpnr'] = 5
taxwcpnrs['Seq363','cpnr'] = 8
taxwcpnrs['Seq36','cpnr'] = 4
taxwcpnrs['Seq44','cpnr'] = 6
taxwcpnrs['Seq82','cpnr'] = 4


#Now that we have finished species level resolution, we move to those that are at genus level
#file from https://rrndb.umms.med.umich.edu/static/download/
allcopynrs = read.csv('rrnDB-5.8_pantaxa_stats_RDP.tsv',sep = '\t')

for (i in which(is.na(taxwcpnrs$cpnr))) {
  
  print(taxwcpnrs[i,'Genus'])
  
  if(!is.na(taxwcpnrs[i,'Genus'])){
    
    
    if(length(which(allcopynrs$name== taxwcpnrs[i,'Genus']))!=0){
      gencpnr = allcopynrs[which(allcopynrs$name == taxwcpnrs[i,"Genus"]),"median"]
      
      taxwcpnrs[i,"cpnr"] = gencpnr
    }
    
  }
}



taxwcpnrs[which(is.na(taxwcpnrs$cpnr)),'cpnr'] = median(taxwcpnrs$cpnr,na.rm = T)



copynrs16s["Clostridium nexile",2] = 4
copynrs16s["Fusobacterium nucleatum",2] = 5
copynrs16s["Bifidobacterium pseudocatenulatum/catenulatum",2] = 5
copynrs16s["Enterococcus faecium/faecalis",2] = 5
copynrs16s["Erysipelatoclostridium ramosum",2] = 5
copynrs16s["Escherichia-Shigella",2] = 7
copynrs16s["Staphylococcus MULSP",2] = 5
copynrs16s["Actinomyces odontolyticus",2] = 3
copynrs16s["Bifidobacterium longum/breve",2] = 5
copynrs16s["Lactobacillus acidophilus/crispatus",2] = 5
copynrs16s["Blautia hansenii",2] = 5
copynrs16s["Eubacterium rectale",2] = 5

