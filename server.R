########Wormtracks server####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###Server

##Required packages
#install.packages("shiny")
#install.packages("shinythemes")


#Load libraries
library(shiny)
library(shinythemes)
library(Biostrings)
library(DT)

shinyServer(function(input, output, session) {
    
    ##########Session functions
    ##Retrieve unique ID for the session
    session_id <- session$token
    
    ##Create temporary folder for unique user
    system(paste("mkdir -p WorkingSpace/users/",session_id,sep=""))
    
    ###On exit, force the remove of directory
    ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
    session$onSessionEnded(function(){
        system(paste("rm -rf WorkingSpace/users/",session_id,sep=""))
    }
    )
    #################################
    
    ###Control panels functions##########################################################################################
    ##Functions needed to generate links between panels
    observeEvent(input$link_to_tabpanel_title, {
        newvalue <- "Designer"
        updateTabsetPanel(session, "panels", newvalue)
    })
    
    observeEvent(input$link_to_tabpanel_downloads, {
        newvalue <- "Downloads"
        updateTabsetPanel(session, "panels", newvalue)
    })
    
    observeEvent(input$link_to_tabpanel_design, {
        newvalue <- "Designer"
        updateTabsetPanel(session, "panels", newvalue)
    })
    
    #########################################################################################################
    ####Database NO selection

    Genes=read.table("WorkingSpace/Gene_DB.tsv",sep="\t", header= FALSE, stringsAsFactors=F)
    MainDB=read.table("WorkingSpace/Main_DB.tsv",sep="\t", header= FALSE, stringsAsFactors=F)
    locusDB=unique(as.character(MainDB[,3]))
    genesDB=unique(as.character(Genes[,1]))
    transDB=as.character(Genes[,2])
    rownames(Genes)=transDB
    
    SizesDB=read.table("WorkingSpace/CDS_sizes.tab",sep="\t",header=F)
    rownames(SizesDB)=as.character(SizesDB[,2])
    
    ##Load data to display for common fluorophores
    fluo=read.table("WorkingSpace/Fluo_sizes.tab", sep="\t", stringsAsFactors=F)
    colnames(fluo)=c("Name","File","Size")
    
    FluoPi=read.table("WorkingSpace/Fluo_master_simplified.tab",sep="\t",stringsAsFactors=F)
    colnames(FluoPi)=c("Plasmid","piRNA","MMG","MMF","Pos","Nuc","GC")
##################################################################################
    
    ##Main search function
    observeEvent(input$actiongenesearch, {
        output$ErrorMessage <- renderText({})
        wbid=""
        mygene = as.character(input$geneinput)

        if(mygene %in% genesDB){
            wbid=mygene
            }
        if(mygene %in% transDB){
            wbid=as.character(unique(Genes[mygene,1]))
            }
        if(mygene %in% locusDB){
            wbid=as.character(MainDB[which(mygene==as.character(MainDB[,3]))[1],1])
            }
        
        if(wbid == ""){
            output$DesignControls <- renderUI({
                HTML("<b>Gene not found</b>")
                })
            output$SelPiTabSummary <- renderUI({})
            output$SelPiTab=renderTable({})
            output$downloadseq <- renderUI({})
            output$SimpleFragment <- renderText({})
            }else{
        output$DesignControls <- renderUI({
            fluidRow(
                
                selectInput("isoform", label = HTML("<b>Isoform
                                                           [<a href=\"\" onclick=\"$('#explain_isoform').toggle(); return false;\">info</a>]
                                                           </b>"), 
                            Genes[which(as.character(Genes[,1])==wbid),2]),
                
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_isoform\">
            WormBase 270 annotations</div></p>
                     "),
                selectInput("selectMM", label = HTML("<b>sg-piRNA specificity
                                                           [<a href=\"\" onclick=\"$('#explain_uniqueness').toggle(); return false;\">info</a>]
                                                           </b>"), 
                            choices = list("Stringent (≥ 4MM)" = 1, "Moderate (≥ 3MM)" = 2, 
                                           "Relaxed (≥ 2MM)" = 3),
                            selected = 1),
                HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness\">
            This option specifies the minimum Hamming distance (i.e., the number of mismatches) to other protein-coding genes.
                                                 </div></p>
                     "),
            #     radioButtons("cluster", label = HTML("Select piRNA cluster
            #                                          [<a href=\"\" onclick=\"$('#explain_cluster').toggle(); return false;\">info</a>]
            #                                          "),
            #                  choices = list("21ur-1224" = 1), selected = 1, width='100%'),
            #     HTML("
            #          <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_cluster\">
            # For the moment, we use the cluster 21ur-1224 as a template to express 6 piRNAis fragments that are antisente to the transcript being targeted
            #                          </div></p>
            #          "),
            radioButtons("Simp_dist", label = HTML("<b>cDNA targeting
                                                               [<a href=\"\" onclick=\"$('#explain_simp_dist').toggle(); return false;\">info</a>]
                                                               </b>"),
            choices = list("Distributed" = 1, "Towards the 3' end" = 2), selected = 1, width='100%', inline= TRUE),
        	HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_simp_dist\">
            Select the location of sg-piRNA binding, i.e., sg-piRNAs distributed uniformly or near the 3' end of the transcript.
                                                 </div></p>
                     "),
            ###Uracil complementary
            radioButtons("Uracil_comp", label = HTML("<b>5' Uracil complementarity
                                                               [<a href=\"\" onclick=\"$('#explain_ura_comp').toggle(); return false;\">info</a>]
                                                               </b>"),
                         choices = list("No preference" = 1, "Cytosine" = 2), selected = 1, width='100%', inline= TRUE),
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_ura_comp\">
            Select the 5' sg-piRNA complementarity. <a href=\"https://www.sciencedirect.com/science/article/pii/S009286741830117X\">Evidence</a> suggests that native piRNAs preferentially match the leading 5' uracil with a cytosine.
                                                 </div></p>
                     "),
            checkboxInput("FlaControl", label = HTML("<b>Negative control
                                                               [<a href=\"\" onclick=\"$('#explain_control').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
            HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_control\">
            Generate negative control piRNAi transgene with sg-piRNAs in the sense orientation (non-silencing).
                                     </div></p>"),
            actionButton("actionPI", label = "Generate piRNAi cluster"),
            
            HTML("<br><br>Please cite Priyadarshini <i>et al.</i>, (2021), Reprogramming the piRNA pathway for multiplexed and transgenerational gene silencing in <i>C. elegans</i>. Nature Methods.")
            
            )
        })
        }
    }, ignoreInit = T)
    
    
    ##Main search function
    observeEvent(input$actionPI, {
        ErrorFlag = 0
        output$ErrorMessage <- renderText({})
        
        matches=as.integer(input$selectMM)
        mm=c(4,3,2)[matches]
        isoform = as.character(input$isoform)
        ControlEx = input$FlaControl
        
        wheretarg=as.integer(input$Simp_dist)
        
        fiveprimecomp=as.integer(input$Uracil_comp)
        
        wbid = as.character(unique(Genes[isoform,1]))
        loc = as.character(unique(Genes[isoform,3]))
        
        file=paste(c("DataBase/",as.character(wbid),"_",as.character(isoform),".txt"),sep="",collapse="")
        
        tab=read.table(file,sep="\t",header=F, stringsAsFactors=F)
        tab[,5]=tab[,4]
        tab[,4]=as.integer((unlist(strsplit(as.character(tab[,2]),";")))[c(FALSE,TRUE)])
        tab[,2]=as.character((unlist(strsplit(as.character(tab[,2]),";")))[c(TRUE,FALSE)])

        Seltab=tab[which(tab[,3]>=mm),]
        
        Seltab = Seltab[which((Seltab[,4]>=30)&(Seltab[,4]<=70)),]
        
        if(fiveprimecomp == 2){
            Seltab = Seltab[which(as.character(Seltab[,5])=="C"),]
            }
        idx=c()
        if(nrow(Seltab)< 6){
            output$ErrorMessage <- renderText({
                paste("Error: Not enough piRNAi fragments were found with the characterisitics described. Try to change to other parameters")
            })
            ErrorFlag=1
        }else{
            Seltab=Seltab[order(Seltab[,1]),]
            
            #Distributed location
            if(wheretarg == 1){
            pos=quantile(Seltab[,1],c(0,.2,.4,.6,.8,1))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[1])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[2])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[3])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[4])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[5])))
            idx=append(idx,which.min(abs(Seltab[,1]-pos[6])))
            }else{
            	#3' targeting
            	idx=append(idx,nrow(Seltab))
            	idx=append(idx,max(which(Seltab[,1] < (Seltab[idx[1],1]-21))))
            	idx=append(idx,max(which(Seltab[,1] < (Seltab[idx[2],1]-21))))
            	idx=append(idx,max(which(Seltab[,1] < (Seltab[idx[3],1]-21))))
            	idx=append(idx,max(which(Seltab[,1] < (Seltab[idx[4],1]-21))))
            	idx=append(idx,max(which(Seltab[,1] < (Seltab[idx[5],1]-21))))
            	
            	}
            
            if(length(which(dist(Seltab[idx,1])<=20)) > 0){
                
                output$ErrorMessage <- renderText({
                    paste("Error: The program selected overlapping piRNAi sites. Try to change the parameters or use Advanced function")
                })
                ErrorFlag=1
                
                }
        }
        
        ##Error for at least 6 piRNAs
        if((length(idx) < 6) | (sum(is.na(Seltab[idx,1]))>0)){
            
            output$ErrorMessage <- renderText({
                paste("Error: Not enough piRNAi sites to create the cluster. Try to change the parameters or use advanced function")
            })
            ErrorFlag=1
            
        }
        
        #Remove if errors()
        if(ErrorFlag == 1){
            output$SelPiTabSummary <- renderUI({ HTML(paste0("<b>Try again!</b>",sep=""))})
            output$SelPiTab=renderTable({})
            output$downloadseq <- renderUI({})
            output$SimpleFragment <- renderText({})
            
            }
        ##Produce outputs
            if(ErrorFlag == 0){
                
                ##Table results
                Pitab=Seltab[idx,c(1,2,4,5)]
                Pistrt=Pitab[,1]
                Piedt=Pitab[,1]+19

                Pitab=cbind(paste(as.integer(Pistrt),"to", as.integer(Piedt)),Pitab[,c(2,3,4)])
                
                colnames(Pitab)=c("Location","Sequence (antisense to target)","%GC","5' complementarity")
                colnames(Pitab)[1]=paste("cDNA location (", as.integer(SizesDB[isoform,3]),"bp long)",sep="")
                Pitab=Pitab[order(Pistrt),]
                
                output$SelPiTabSummary <- renderUI({ HTML(paste0("<b>Synthetic piRNAs</b>",sep=""))})
                output$SelPiTab=renderTable(Pitab)
                
                ##Ape output
                output$downloadseq <- renderUI({
                    ##If control experiment, invert sequences
                    if(ControlEx){
                        Seltab[idx[1],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[1],2]))))
                        Seltab[idx[2],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[2],2]))))
                        Seltab[idx[3],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[3],2]))))
                        Seltab[idx[4],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[4],2]))))
                        Seltab[idx[5],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[5],2]))))
                        Seltab[idx[6],2]=as.character(reverseComplement(DNAString(as.character(Seltab[idx[6],2]))))
                        }
                    
                    uno="cgcgcttgacgcgctagtcaactaacataaaaaaggtgaaacattgcgaggatacatagaaaaaacaatacttcgaattcatttttcaattacaaatcctgaaatgtttcactgtgttcctataagaaaacattgaaacaaaatattAagT"
                    uno=tolower(uno)
                    seq1=as.character(Seltab[idx[1],2])
                    dos="ctaattttgattttgattttgaaatcgaatttgcaaatccaattaaaaatcattttctgataattagacagttccttatcgttaattttattatatctatcgagttagaaattgcaacgaagataatgtcttccaaatactgaaaatttgaaaatatgtt"
                    dos=tolower(dos)
                    seq2=as.character(reverseComplement(DNAString(as.character(Seltab[idx[2],2]))))
                    tres="AttGccagaactcaaaatatgaaatttttatagttttgttgaaacagtaagaaaatcttgtaattactgtaaactgtttgctttttttaaagtcaacctacttcaaatctacttcaaaaattataatgtttcaaattacataactgtgt"
                    tres=tolower(tres)
                    seq3= as.character(reverseComplement(DNAString(as.character(Seltab[idx[3],2]))))
                    cuatro="ActgtagagcttcaatgttgataagatttattaacacagtgaaacaggtaatagttgtttgttgcaaaatcggaaatctctacatttcatatggtttttaattacaggtttgttttataaaataattgtgtgatggatattattttcagacctcatactaatctgcaaaccttcaaacaatatgtgaagtctactctgtttcactcaaccattcatttcaatttggaaaaaaatcaaagaaatgttgaaaaattttcctgtttcaacattatgacaaaaatgttatgattttaataaaaaCaaT"
                    cuatro=tolower(cuatro)
                    seq4=as.character(Seltab[idx[4],2])
                    cinco="ttctgtttttcttagaagtgttttccggaaacgcgtaattggttttatcacaaatcgaaaacaaacaaaaatttttttaattatttctttgctagttttgtagttgaaaattcactataatcatgaataagtgagctgcccaagtaaacaaagaaaatttggcagcggccgacaactaccgggttgcccgatttatcagtggagga"
                    cinco=tolower(cinco)
                    seq5= as.character(reverseComplement(DNAString( as.character(Seltab[idx[5],2]))))
                    seis="AtcTaatgtgatgtacacggttttcatttaaaaacaaattgaaacagaaatgactacattttcaaattgtctatttttgctgtgtttattttgccaccaacaaT"
                    seis=tolower(seis)
                    seq6=as.character(Seltab[idx[6],2])
                    siete="tcaatctagtaaactcacttaatgcaattcctccagccacatatgtaaacgttgtatacatgcagaaaacggttttttggttttaatgggaacttttgacaaattgttcgaaaatcttaagctgtcccatttcagttgggtgatcgattt"
                    siete=tolower(siete)

                    xtracom=paste("Standard piRNA transgene (ClusterE). Parameters: Gene = ",wbid,"; Isoform = ",isoform,"; At least ",mm," mismatches")
                    binrev=c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
                    if(ControlEx){
                       xtracom = append(xtracom,"Control experiment: piRNAi Sequences have been reverse complemented. THIS FRAGMENT WONT SILENCE THE SELECTED GENE")
                       binrev=!binrev
                    }

                    Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                    
                    toadd="WorkingSpace/Piconst2.txt"
                    
                    pats= c(seq1,seq2,seq3,seq4,seq5,seq6)
                    fwdc= c(rep("#00ff00",6))
                    revc= c(rep("#ff0000",6))
                    tooltis= paste("piRNA",1:6)

                    writeLines(PasteApe(paste(wbid,"_Standard_ClusterE_",sep="",collapse=""),Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"Caenorhabditis"),paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""))
                    
                    output$SimpleFragment <- renderText({
                        paste("Recoded standard piRNA transgene (ClusterE)\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                    })
                    
                    downloadButton('DownApeOut', 'Download annotated genbank file')
                })
                
                
            }
        
    }, ignoreInit = T)
    
    ###Advanced construct######
    ####Dirty coded addition of extra clusters
    observeEvent(input$actionconstruct, {
        AdvancedErrorFlag=0
        output$AdvancedErrorMessage <- renderText({})
        clust=as.character(input$clustercon)

        if(clust == "1"){

        pipi1=as.character(input$piRNAseq1_1)
        pipi2=as.character(input$piRNAseq2_1)
        pipi3=as.character(input$piRNAseq3_1)
        pipi4=as.character(input$piRNAseq4_1)
        pipi5=as.character(input$piRNAseq5_1)
        pipi6=as.character(input$piRNAseq6_1)
        
        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 6 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="cgcgcttgacgcgctagtcaactaacataaaaaaggtgaaacattgcgaggatacatagaaaaaacaatacttcgaattcatttttcaattacaaatcctgaaatgtttcactgtgttcctataagaaaacattgaaacaaaatattAagT"
                uno=tolower(uno)
                seq1=as.character(pipi1)
                dos="ctaattttgattttgattttgaaatcgaatttgcaaatccaattaaaaatcattttctgataattagacagttccttatcgttaattttattatatctatcgagttagaaattgcaacgaagataatgtcttccaaatactgaaaatttgaaaatatgtt"
                dos=tolower(dos)
                seq2=as.character(reverseComplement(DNAString(as.character(pipi2))))
                tres="AttGccagaactcaaaatatgaaatttttatagttttgttgaaacagtaagaaaatcttgtaattactgtaaactgtttgctttttttaaagtcaacctacttcaaatctacttcaaaaattataatgtttcaaattacataactgtgt"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="ActgtagagcttcaatgttgataagatttattaacacagtgaaacaggtaatagttgtttgttgcaaaatcggaaatctctacatttcatatggtttttaattacaggtttgttttataaaataattgtgtgatggatattattttcagacctcatactaatctgcaaaccttcaaacaatatgtgaagtctactctgtttcactcaaccattcatttcaatttggaaaaaaatcaaagaaatgttgaaaaattttcctgtttcaacattatgacaaaaatgttatgattttaataaaaaCaaT"
                cuatro=tolower(cuatro)
                seq4=as.character(pipi4)
                cinco="ttctgtttttcttagaagtgttttccggaaacgcgtaattggttttatcacaaatcgaaaacaaacaaaaatttttttaattatttctttgctagttttgtagttgaaaattcactataatcatgaataagtgagctgcccaagtaaacaaagaaaatttggcagcggccgacaactaccgggttgcccgatttatcagtggagga"
                cinco=tolower(cinco)
                seq5= as.character(reverseComplement(DNAString( as.character(pipi5))))
                seis="AtcTaatgtgatgtacacggttttcatttaaaaacaaattgaaacagaaatgactacattttcaaattgtctatttttgctgtgtttattttgccaccaacaaT"
                seis=tolower(seis)
                seq6=as.character(pipi6)
                siete="tcaatctagtaaactcacttaatgcaattcctccagccacatatgtaaacgttgtatacatgcagaaaacggttttttggttttaatgggaacttttgacaaattgttcgaaaatcttaagctgtcccatttcagttgggtgatcgattt"
                siete=tolower(siete)
                
                xtracom=paste("Standard piRNA transgene (ClusterE) via Advanced Search.")
                binrev=c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                
                toadd="WorkingSpace/Piconst2.txt"
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6)
                fwdc= c(rep("#00ff00",6))
                revc= c(rep("#ff0000",6))
                tooltis= paste("piRNA",1:6)
                
                writeLines(PasteApe("Standard_ClusterE_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"Caenorhabditis"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded standard piRNA transgene (ClusterE)\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        ##Second cluster
        
        if(clust == "2"){

        pipi1=as.character(input$piRNAseq1_2)
        pipi2=as.character(input$piRNAseq2_2)
        pipi3=as.character(input$piRNAseq3_2)
        pipi4=as.character(input$piRNAseq4_2)
        pipi5=as.character(input$piRNAseq5_2)
        pipi6=as.character(input$piRNAseq6_2)
        
        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 6 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="ttcgtggtgcacttatctttctccttcaaattgaaaactcagtttttaattatagtcaaatctcttttgctgacaggtccaaagtactttattatttcatattatataaaattcattctcgaatttatttataaattttcgctgagtcaa"
                uno=tolower(uno)
                seq1=as.character(reverseComplement(DNAString(as.character(pipi1))))
                dos="ActgactaacaaaaacccctgtcaatttacttgtaatgtgaaactgtatcggtttcatattatctatgattcgagtacattgtttcaaattcaaT"
                dos=tolower(dos)
                seq2=as.character(pipi2)
                tres="attttacgctggtttgaaaatttgaaatattccaaaataaattatttagttttcgttttttgtacattgtcataaaacattttggttttttttaaca"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="Acttaaaatcaaaaattgttacactttataacagttcattgaaactgaaaattattttcttttcccaaataataataccatcaaatgtcgtggtgtactcatcttttccttttcttcttttttttcaatttctccttcaaatctctacacactcttcactgccaatctttttttctttccttatccaaaagcacacttttgtgcagagtaaataatgcactttgtgaaaaaaaaactatttttaaaactgtatttttttaagtttggcaatttttgagaaaatttcaacaaaatctgatatagattggaatttaaatggttcaaatttg"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="AtctattcaaagttttattcgaagtttttaacagacacttgaaacagtgtaataattttctgacaaaaattaaaacaaatgttactactttgcttttcttactttatccgttttttatcacccttatttttcagtcaaccctagcaacgttaccgacggaatcggtaggactacaccgactgcatcaaatttgggaagaagccgtgagaatttgagtttcaatcaacaccgcccagaccatccttcgtcatattttgatagtttggagcatggtgagcattttatattaaaacagttgttttggtgttcatattactaatgtctgaatactaacttgcattaaaattggaaattaaaaaaattactgtttctcaaaagtattttcaatacctatatttttttgctacagT"
                cinco=tolower(cinco)
                seq5= as.character(pipi5)
                seis="caatattttcaaatattttataccagatttttcgaaaaagttgaattttcaattaacaataacgcatttatgcatttttcactcttttttgagatttaatgctgaaaaaatagttctgaaaatgacaaaagttatgttttcaatattttttatcaaactaaatttatttaatttgttaactgttgcttttttgtttttcttcaagt"
                seis=tolower(seis)
                seq6=as.character(reverseComplement(DNAString(as.character(pipi6))))
                siete="Atcttcgaagcaacttatttgatgttttataaacgacctgaaacatactggtgatgcccaataatgttttttttaaatttagtctcgtgaaaaaaataaaattaaaacagaaaattacatttgcgccgaagaaacttaagatctggaactt"
                siete=tolower(siete)
                
                xtracom=paste("Recoded ClusterG via Advanced Search.")
                binrev=c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6)
                fwdc= c(rep("#00ff00",6))
                revc= c(rep("#ff0000",6))
                tooltis= paste("piRNA",1:6)
                
                writeLines(PasteApe("ClusterG_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"Caenorhabditis"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded piRNA ClusterG\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        ###Third cluster

        if(clust == "3"){

        pipi1=as.character(input$piRNAseq1_3)
        pipi2=as.character(input$piRNAseq2_3)
        pipi3=as.character(input$piRNAseq3_3)
        pipi4=as.character(input$piRNAseq4_3)
        pipi5=as.character(input$piRNAseq5_3)
        pipi6=as.character(input$piRNAseq6_3)
        pipi7=as.character(input$piRNAseq7_3)
        
        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)| (nchar(pipi7) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 7 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)| (nchar(pipi7) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="agggtggtcgcatagaagttgggcgcacttcatttacaaaaatatgttcaaattttgtgatttcatgttcaggcattttgttttgatgataaaacatagtgtgactgtttttatatgtttataaaatgtcttatgattaaacaagatcaaT"
                uno=tolower(uno)
                seq1=as.character(pipi1)
                dos="aaatgtcggtttatttcacactgataggaatttttcaaaaaaatatatcagcaaagtactgtattaaaatgtgaaaatctcataaaaagtttaagtttca"
                dos=tolower(dos)
                seq2=as.character(reverseComplement(DNAString(as.character(pipi2))))
                tres="ATtgtggcaaatagattttgacaatttttatcaaattcatgaaacagtagaattttttccaaaaaactcacaaaataaatacgaatttcaatttgcccactttatcaaataaatgtttacacaaaagtaggccgtgcaacgcgcctatcctagatgctacattccttggttttgagttgtgaaacgttggaataatgtacttcattttgtgacttacttttttgtaaggcaattgttttttttatttaataaaagtactttcctaaattcaaatatcaaatttgtcttcatttttgtgacaagtaaaa"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="Agcttttttagaaaaaaaaagtcagtatttaagaacaatttgaaacagttcaatttttcagtgtaatttccaaccagaactttttgagtaaataattacagaaaacttatttagaaaataggactaaataatgcaaatattttccggactggcatttataaataagcaagtataag"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="AttgaggtaattttaaaaagcatacatataagcaagtcgtgaaacagtcgtttaaattttatttttcaaaaagttataacgcgacagcagtttcatctgtttcatattccctatttgttgaaatttgagacgtattttacgaT"
                cinco=tolower(cinco)
                seq5= as.character(pipi5)
                seis="tgccatccgaatcttgaactttgtatcaattgttcacatttttttccaaaaacgtattaactcactttca"
                seis=tolower(seis)
                seq6=as.character(reverseComplement(DNAString(as.character(pipi6))))
                siete="AtcttttgtttttaacaaagagatcatatatactcattgagaaacagtacattttttgaaagtacatttgccttggtcaaatatataaagttgacaaaagtttaaaaatgtttccaaaagttaattaaaaaatcaaatttattctagctcaaT"
                siete=tolower(siete)
                seq7=as.character(pipi7)
                ocho="cgtcaagtgatcaaaccatcatttttttcagattaagacctgatttgtcagtgaattgaaaaaaacgtgttcattgcgtgtttcgcattttttatatataaaaaagcaagtttcggcggcaataacgaagtattcccaacagatcaatag"
                ocho=tolower(ocho)
                
                xtracom=paste("Recoded ClusterO via Advanced Search.")
                binrev=c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete,seq7,ocho),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6,seq7)
                fwdc= c(rep("#00ff00",7))
                revc= c(rep("#ff0000",7))
                tooltis= paste("piRNA",1:7)
                
                writeLines(PasteApe("ClusterO_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"Caenorhabditis"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded piRNA ClusterO\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        ###Fourth
        
        if(clust == "4"){

        pipi1=as.character(input$piRNAseq1_4)
        pipi2=as.character(input$piRNAseq2_4)
        pipi3=as.character(input$piRNAseq3_4)
        pipi4=as.character(input$piRNAseq4_4)
        pipi5=as.character(input$piRNAseq5_4)
        pipi6=as.character(input$piRNAseq6_4)
        pipi7=as.character(input$piRNAseq7_4)
        pipi8=as.character(input$piRNAseq8_4)

        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)| (nchar(pipi7) == 0)| (nchar(pipi8) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 8 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)| (nchar(pipi7) != 20)| (nchar(pipi8) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,pipi8,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,pipi8),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="caaaaaacaatacgtcccttatcttctggaatcagctcattgtgctcatcggagctatccgcaccgtcaactatactcgctagatcttccgtgttctgatcttgagtgtatagtggaggggggtcaacctgaaatttcagatttttgttg"
                uno=tolower(uno)
                seq1=as.character(reverseComplement(DNAString(as.character(pipi1))))
                dos="ActgttttagaagtgatgagtcttattataataacttgttgaaactgtggatttatattttttaaaaattaccggcgaaattgattcataatctcttattaccatagttaaagtctctagaataagcacaaaactactaaagtttgtaaaataattgaatatgccacaactgataagagactttttcctcttatcagcataaagtccaaagcgataaaattcaaaagagacaagtacaaatgtatattaatctgctttgttggaaaaaaattaaactttttatctaaacctgtcattgatccaaaagattaagtttcctgcaaaattgtttcgaaatattattgtgattgaaacttttgactttttcaacttatcaataagtcattggcttaagataaagtaatcaaT"
                dos=tolower(dos)
                seq2=as.character(pipi2)
                tres="cgcgctcagcactcaatttctgcccaaatagtt"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="Attgagtatctaaatgaaaacctaaaatatgaacagttagaaacaggaaatttttgaaaagttaaaaaacaacctatacaattaatttccaagaaaaatttaacaatcgattttcatttctgaaatcccaaaatcggtgaattcttgatgaaaatgcatttgaaaatacaattttgtttta"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="Atctaattagatatgcaagcctaatatttgtatcattcttgaaactgtaaataaaaaatgtttgcaaaaaaaatcaattttttagcgaatgttaacataaaaccttaaatttttctgggttttgaccgtttctcatatttcaaaa"
                cinco=tolower(cinco)
                seq5= as.character(reverseComplement(DNAString(as.character(pipi5))))
                seis="AtcgaT"
                seis=tolower(seis)
                seq6=as.character(pipi6)
                siete="aaaaaatatgctgaaacgtgattgctttttgtgcttttttatacaagtttgcaatcgcacaaatcatatgaaaaattattaagcacgcttaaactatgtgatctgaaatacgaaaactagtatacgttaaacaggaaaaaaaaatcaactgtttcaaaatttgtgtttaatcaaattaaagttgctattccgaT"
                siete=tolower(siete)
                seq7=as.character(pipi7)
                ocho="taatagttaattttcaaaatagaaagttttaaaacatcctgtttcggtatgctgatttttacagactccactttgtagttaacT"
                ocho=tolower(ocho)
                seq8=as.character(pipi8)
                nueve="ggacagaaatattgatattttgccagttaccaggaaaataaattattctttgcaacatctgactttaagaataaaaactcacaaattccttttccatttctgaaatattttagtgtcctcttcccgcaaccactccctgtaaatcgaaaa"
                nueve=tolower(nueve)
                
                xtracom=paste("Recoded ClusterF via Advanced Search.")
                binrev=c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete,seq7,ocho,seq8,nueve),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8)
                fwdc= c(rep("#00ff00",8))
                revc= c(rep("#ff0000",8))
                tooltis= paste("piRNA",1:8)
                
                writeLines(PasteApe("ClusterF_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"Caenorhabditis"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded piRNA ClusterF\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

         ###Fifth
        
        if(clust == "5"){

        pipi1=as.character(input$piRNAseq1_5)
        pipi2=as.character(input$piRNAseq2_5)
        pipi3=as.character(input$piRNAseq3_5)
        pipi4=as.character(input$piRNAseq4_5)
        pipi5=as.character(input$piRNAseq5_5)
        pipi6=as.character(input$piRNAseq6_5)
        pipi7=as.character(input$piRNAseq7_5)
        pipi8=as.character(input$piRNAseq8_5)

        ##Check for size
        if((nchar(pipi1) == 0) | (nchar(pipi2) == 0) | (nchar(pipi3) == 0) | (nchar(pipi4) == 0) | (nchar(pipi5) == 0) | (nchar(pipi6) == 0)| (nchar(pipi7) == 0)| (nchar(pipi8) == 0)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Please pick at least 8 synthetic piRNAs")
            })
        }
        
        if(AdvancedErrorFlag == 0){
        if((nchar(pipi1) != 20) | (nchar(pipi2) != 20) | (nchar(pipi3) != 20) | (nchar(pipi4) != 20) | (nchar(pipi5) != 20) | (nchar(pipi6) != 20)| (nchar(pipi7) != 20)| (nchar(pipi8) != 20)){
            AdvancedErrorFlag=1
            output$AdvancedErrorMessage <- renderText({
                paste("Error: piRNAi sequences should be 20bp long")
            })
            }}
        
        #Check for input characters
        toto=paste(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,pipi8,sep="",collapse="")
        if((AdvancedErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(toto))) != 0)){ ##Check for strange non ATCG characters
            output$AdvancedErrorMessage <- renderText({
                paste("Error: Unrecognized characters in piRNAi sequences")
            })
            AdvancedErrorFlag=1
        }
        
        #Checkfor GC
        if(AdvancedErrorFlag == 0){
        Gcvals=sapply(c(pipi1,pipi2,pipi3,pipi4,pipi5,pipi6,pipi7,pipi8),CalculateGC)
        if((sum(Gcvals <.3)+sum(Gcvals >.7))>0){
            output$AdvancedErrorMessage <- renderText({
                paste("Warning: some sequences have a GC content lower to 30% or higher to 70%")
            })
        }
        }
        
        ##Main Routine
        if(AdvancedErrorFlag == 0){
            output$downloadconstruct <- renderUI({
            
                uno="gaatttgttattttctatcatattgacaaaacaaaaaaacatattcaaaa"
                uno=tolower(uno)
                seq1=as.character(reverseComplement(DNAString(as.character(pipi1))))
                dos="aTtaccctaaattttaaaaaaataattatacaaacgtagtgaaacagcagattaattttccattaatatacaacaaaagtttgtttttaattgaaatttcagattttaagtaaat"
                dos=tolower(dos)
                seq2=as.character(reverseComplement(DNAString(as.character(pipi2))))
                tres="aTcactaaagaattttcaagtaaactatataaacgtagtaaaacagcagaaatatttttagcattttt"
                tres=tolower(tres)
                seq3= as.character(reverseComplement(DNAString(as.character(pipi3))))
                cuatro="atcaaagcacatatttttgacggaaactatatgaacttggtgaaacagtagaaattggaaaaaaaaatagtgtttcacaacttttatgaaagttttgattttgttttctaaattaaaaatatctttgagctatgataaaaattgatacaatcacttgtgtttcaactgcaatcatgttaattacggaatccgtgagctatttgaaaattgaaattatt"
                cuatro=tolower(cuatro)
                seq4=as.character(reverseComplement(DNAString(as.character(pipi4))))
                cinco="attgaagtacaaaattcaacataaaattatatcaacaagttgaaacagttcttacattttttgaaaattgtcgactttttttttacattacccatgtaatttttaaatcaaatttaaaattattcgtattaccttagtaatgggattatttttgtttacgtgattgtcagttgaaaaattgattttttaatggtgacagggatctgtttcgttaagttactaaaataagaaccaattcctccaacttggAt"
                cinco=tolower(cinco)
                seq5= as.character(pipi5)
                seis="aataaattattagtgattgtaagtaaaattatttcaaaatttttctacatgtatacactacatgttcacgagcataaagttaaaaataga"
                seis=tolower(seis)
                seq6=as.character(reverseComplement(DNAString(as.character(pipi6))))
                siete="atcaagttgctaaatttttaatggaattaaatgaagatggtgaaactgtagaaaatttttaaaaatttttttgctgtttcacaattaaaaaaaataattcaggcaaaatttcaattgAt"
                siete=tolower(siete)
                seq7=as.character(pipi7)
                ocho="gtattactatttttatcgactcaatcaccgatttagaaatttaaattcaagttttcaaactatgtataattttgtttaaaggtctttttgaattttttaacttgcaataaaagctgccatggaaattatcttttagcaatatgcacacttttaacaaggtaggcaggcgttttcgtacctacactgcagatcatatatattaattccgccaataaaagatcttgaaattaaaaaaaaactactaatttcaagacataattgaaaaatttagattat"
                ocho=tolower(ocho)
                seq8=as.character(reverseComplement(DNAString(as.character(pipi8))))
                nueve="atcagaaaacacttttttttaacgttcttatgtgaacatttagaaacagtg"
                nueve=tolower(nueve)
                
                xtracom=paste("Recoded piRNA cluster centered on 21ur-5764 via Advanced Search.")
                binrev=c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
                
                Compseq=paste(c(uno,seq1,dos,seq2,tres,seq3,cuatro,seq4,cinco,seq5,seis,seq6,siete,seq7,ocho,seq8,nueve),sep="",collapse="")
                
                toadd=""
                
                pats= c(seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8)
                fwdc= c(rep("#00ff00",8))
                revc= c(rep("#ff0000",8))
                tooltis= paste("piRNA",1:8)
                
                writeLines(PasteApe("21ur-5764_",Compseq,pats,fwdc,revc,tooltis,xtracom,toadd,binrev,"Caenorhabditis"),paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""))
                
                output$AdvancedFragment <- renderText({
                    paste("Recoded 21ur-5764 piRNA cluster\n",paste(Compseq,sep="",collapse=""),sep="",collapse="")
                })
                
            downloadButton('DownConOut', 'Download annotated genbank file')
            
                })
            
            }
        }

        }, ignoreInit = T)
    
    ###Advanced searchform
    ###Only upon retrieval of a gene in database new form will appear
    observeEvent(input$actionAdvsearch, {
        wbid=""
        mygene = as.character(input$Advancedgeneinput)
        
        if(mygene %in% genesDB){
            wbid=mygene
        }
        
        if(mygene %in% transDB){
            wbid=as.character(unique(Genes[mygene,1]))
        }
        
        if(mygene %in% locusDB){
            wbid=as.character(MainDB[which(mygene==as.character(MainDB[,3]))[1],1])
        }
        
        if(wbid == ""){
            output$AdvDesignControls <- renderUI({
                HTML("<b>Gene not found</b>")
            })
            
            output$AdvDesignControls <- renderUI({})
            
        }else{
            
            ##Control for table
            output$AdvDesignControls <- renderUI({
                fluidRow(
                    
                    column(width = 2,selectInput("AdvIsoform", label = HTML("<b>Isoform
                                                           [<a href=\"\" onclick=\"$('#explain_isoform_advanced').toggle(); return false;\">info</a>]
                                                           </b>"), 
                                                 Genes[which(as.character(Genes[,1])==wbid),2])),
                    
                    column(width = 3, selectInput("AdvSelectMM", label = HTML("<b>sg-piRNA specificity
                    [<a href=\"\" onclick=\"$('#explain_uniqueness_advanced').toggle(); return false;\">info</a>]
                                                           </b>"),
                                                  choices = list("Stringent (≥ 4MM)" = 1, "Moderate (≥ 3MM)" = 2, 
                                                                 "Relaxed (≥ 2MM)" = 3),
                                                                      selected = 1)),
                    column(width = 3, sliderInput("Posslider", label = HTML("<b>Relative position in cDNA (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_Posgene').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                            "),
                                                  0, 100, c(0, 100), step = 5)),
                    column(width = 2, sliderInput("Gcslider", label = HTML("<b>GC content (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_GCcont').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                           ")
                                                  ,0, 100, c(30, 70), step = 5)),
                    column(width = 2, checkboxGroupInput("CompGroup", label = HTML("<b>5' Uracil complementarity (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_5p_Comp').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                           "),
                                                  choices=c("A","T","C","G"),
                                                  selected=c("A","T","C","G"),
                                                  inline=T)),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_isoform_advanced\">
            WormBase 270 annotations</div></p>
                     "),

                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness_advanced\">
            This option specifies the minimum Hamming distance (i.e., the number of mismatches) to other protein-coding genes.
                                                             </div></p>
                     "),
                    
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_Posgene\">
            Select the location of sg-piRNA binding, i.e., sg-piRNAs distributed uniformly or near the 3' end of the transcript. 
                                                 </div></p>
                     "),
                    
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_GCcont\">
            We recommend piRNAi fragments with GC content between 30 to 70%.
                                                 </div></p>
                     "),
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_5p_Comp\">
            Select the 5' sg-piRNA complementarity. <a href=\"https://www.sciencedirect.com/science/article/pii/S009286741830117X\">Evidence</a> suggests that native piRNAs preferentially match the leading 5' uracil with a cytosine.
                                                 </div></p>
                     ")
                    )
            })
            
            }
        
        
        }, ignoreInit = T)
    
    
    ##Advance search but for common fluorophores
    observeEvent(input$ActionFluoSearch, {
            ##Control for table
            output$AdvDesignControls <- renderUI({
                fluidRow(
                    column(width = 3, selectInput("FluoSelectMMG", label = HTML("<b>sg-piRNA genomic specificity
                    [<a href=\"\" onclick=\"$('#explain_uniqueness_fluo_g').toggle(); return false;\">info</a>]
                                                           </b>"),
                    choices = list("Stringent (≥ 4MM)" = 1, "Moderate (≥ 3MM)" = 2, 
                                   "Relaxed (≥ 2MM)" = 3),
                    selected = 1)),
                    
                    column(width = 2, selectInput("FluoSelectMMF", label = HTML("<b>sg-piRNA fluorophore specificity
                    [<a href=\"\" onclick=\"$('#explain_uniqueness_fluo_f').toggle(); return false;\">info</a>]
                                                           </b>"),
                    choices = list("No filter" = 1,"Stringent (≥ 4MM)" = 2, "Moderate (≥ 3MM)" = 3, 
                                   "Relaxed (≥ 2MM)" = 4),
                    selected = 1)),
                    column(width = 3, sliderInput("FluoPosslider", label = HTML("<b>Relative position in cDNA (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_Posgene_fluo').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                            "),
                                                  0, 100, c(0, 100), step = 5)),
                    column(width = 2, sliderInput("FluoGcslider", label = HTML("<b>GC content (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_GCcont_fluo').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                           ")
                                                  ,0, 100, c(30, 70), step = 5)),
                    column(width = 2, checkboxGroupInput("FluoCompGroup", label = HTML("<b>5' Uracil complementarity (%)
                                                                            [<a href=\"\" onclick=\"$('#explain_5p_Comp_fluo').toggle(); return false;\">info</a>]
                                                                            </b>
                                                                           "),
                                                         choices=c("A","T","C","G"),
                                                         selected=c("A","T","C","G"),
                                                         inline=T)),
                    br(),
                    br(),
                    br(),
                    br(),
                    br(),
                    HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness_fluo_g\">
            This option specifies the minimum Hamming distance (i.e., the number of mismatches) to other protein-coding genes.
            </div></p>
                     "),
            
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_uniqueness_fluo_f\">
            This option specifies the minimum Hamming distance (i.e., the number of mismatches) to common fluorescent proteins.
                                                             </div></p>
                     "),
            
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_Posgene_fluo\">
            Select the location of sg-piRNA binding, i.e., sg-piRNAs distributed uniformly or near the 3' end of the transcript. 
                                                 </div></p>
                     "),
            
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_GCcont_fluo\">
            We recommend piRNAi fragments with GC content between 30 to 70%.
                                                 </div></p>
                     "),
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_5p_Comp_fluo\">
            Select the preferential 5' piRNAi complementarity. <a href=\"https://www.sciencedirect.com/science/article/pii/S009286741830117X\">Evidence</a> suggest that native piRNAs binds preferentially to sites ending in Cytosine.
                                                 </div></p>
                     ")
                )
            })

    }, ignoreInit = F)
    
    ##Observe changes to fluo parameters
    observeEvent({
        ##Parameters of previous options
        input$FluoSelectMMG
        input$FluoSelectMMF
        input$FluoPosslider
        input$FluoGcslider
        input$FluoCompGroup
        input$ActionFluoSearch
        }, {
   
            ##Initial table         
            myfluo =  as.character(fluo[which(as.character(input$FluoInput) == fluo$Name),2])
            lengcdna = as.numeric(fluo[which(as.character(input$FluoInput) == fluo$Name),3])
            
            datatab = FluoPi[which(as.character(FluoPi$Plasmid) == myfluo),]
            
            matches=as.integer(input$FluoSelectMMG)
            mmg=c(4,3,2)[matches]
            matches=as.integer(input$FluoSelectMMF)
            mmf=c(0,4,3,2)[matches]
            
            minGC=input$FluoGcslider[1]
            maxGC=input$FluoGcslider[2]
            minPos=input$FluoPosslider[1]/100
            maxPos=input$FluoPosslider[2]/100
            
            ##5prime comp
            leadingnuc=as.character(input$FluoCompGroup)
            
            datatab = datatab[which(as.character(datatab$Nuc) %in% leadingnuc),]
            
            datatab = datatab[which(datatab$MMG>=mmg),]
            datatab = datatab[which(datatab$MMF>=mmf),]
            
            datatab = datatab[which((datatab$GC>=minGC)&(datatab$GC<=maxGC)),]
            
            CDSlong=as.integer(lengcdna)
            relpos=as.numeric(datatab$Pos)
            
            relpos=relpos/CDSlong
            
            datatab = datatab[which((relpos>=minPos)&(relpos<=maxPos)),]

            
                output$AllPiTab <- DT::renderDataTable({
                    datatab=datatab[order(datatab$Pos),]
                    
                    Pdata=data.frame(
                        Location = paste(as.integer(datatab$Pos), "to", as.integer(datatab$Pos)+19),
                        Sequence = datatab[,2],
                        GCcontent = datatab$GC,
                        Complementarity = datatab$Nuc,
                        Select = shinyInput(actionButton, as.character(datatab[,2]), 'button_', label = "Add to contruct", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
                        stringsAsFactors = FALSE,
                        row.names = 1:nrow(datatab)
                    )
                    
                    colnames(Pdata)[1]=paste("cDNA location ","(",(lengcdna),"bp long)",sep="")
                    colnames(Pdata)[2]="Sequence (antisense to target)"
                    colnames(Pdata)[4]="5' Complementarity"
                    rownames(Pdata)=1:nrow(Pdata)
                    Pdata
                    #},server = FALSE, escape = FALSE, selection = 'none'))
                },server = FALSE, escape = FALSE, selection = 'none')
            
            
            
        
        },ignoreInit = F)
    
    ###Observe function for adding fragments based on piRNA table
    ##############################################################
    
    observeEvent(
        {
        #Track any change to parameters
        input$AdvIsoform
        input$Gcslider
        input$Posslider
        input$AdvSelectMM
        input$CompGroup
        }
        ,{
        ##Now design table; NOt sure if it will work as table should be most of the time be out of observe functions
        ADVisoform = as.character(input$AdvIsoform)
        
        wbid = as.character(unique(Genes[ADVisoform,1]))
        loc = as.character(unique(Genes[ADVisoform,3]))

        file=paste(c("DataBase/",as.character(wbid),"_",as.character(ADVisoform),".txt"),sep="",collapse="")
        
        tab=read.table(file,sep="\t",header=F, stringsAsFactors=F)
        tab[,5]=tab[,4]
        tab[,4]=as.integer((unlist(strsplit(as.character(tab[,2]),";")))[c(FALSE,TRUE)])
        tab[,2]=as.character((unlist(strsplit(as.character(tab[,2]),";")))[c(TRUE,FALSE)])
        
        ##Partial patch to solve for when nopiRNA exist within the input ranges
        datatab = tab
        matches=as.integer(input$AdvSelectMM)
        mm=c(4,3,2)[matches]
        minGC=input$Gcslider[1]
        maxGC=input$Gcslider[2]
        minPos=input$Posslider[1]/100
        maxPos=input$Posslider[2]/100
        
        ##5prime comp
        leadingnuc=as.character(input$CompGroup)
        
        datatab = datatab[which(as.character(datatab[,5]) %in% leadingnuc),]
        
        datatab = datatab[which(datatab[,3]>=mm),]
        
        datatab = datatab[which((datatab[,4]>=minGC)&(datatab[,4]<=maxGC)),]

        CDSlong=as.integer(SizesDB[ADVisoform,3])
        relpos=datatab[,1]

        relpos=relpos/CDSlong

        datatab = datatab[which((relpos>=minPos)&(relpos<=maxPos)),]
        
        if( (nrow(datatab) == 0) | (ncol(datatab) == 0) ){
            output$AllPiTab <- DT::renderDataTable({ data.frame(Error=c("There is no piRNAi target sites with the given parameters")) })
            }else{   
        #output$AllPiTab <- DT::renderDataTable(DT::datatable({
        output$AllPiTab <- DT::renderDataTable({

            lengcdna=CDSlong
            datatab=datatab[order(datatab[,1]),]

            Pdata=data.frame(
                Location = paste(as.integer(datatab[,1]), "to", as.integer(datatab[,1])+19),
                Sequence = datatab[,2],
                GCcontent = datatab[,4],
                Complementarity = datatab[,5],
                Select = shinyInput(actionButton, as.character(datatab[,2]), 'button_', label = "Add to contruct", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
                stringsAsFactors = FALSE,
                row.names = 1:nrow(datatab)
                )
            
            colnames(Pdata)[1]=paste("cDNA location ","(",(lengcdna),"bp long)",sep="")
            colnames(Pdata)[2]="Sequence (antisense to target)"
            colnames(Pdata)[4]="5' Complementarity"
            rownames(Pdata)=1:nrow(Pdata)
            Pdata
        #},server = FALSE, escape = FALSE, selection = 'none'))
        },server = FALSE, escape = FALSE, selection = 'none')
    }
        }, ignoreInit = F)
    
    ##Handle shiny to add dynamic button
    shinyInput <- function(FUN, seqs, id, ...) {
        inputs <- character(length(seqs))
        for (i in 1:length(seqs)) {
            inputs[i] <- as.character(FUN(paste0(id, seqs[i]), ...))
        }
        inputs
    }
    
    #Handle id-seq of dynamic button
    observeEvent(input$select_button, {
        clust=as.character(input$clustercon)
        fill=0
        selectedSeq <- as.character(strsplit(input$select_button, "_")[[1]][2])
        
        if(clust == "1"){
        if(fill == 0){
        if(as.character(input$piRNAseq1_1)==""){
            updateTextAreaInput(session, "piRNAseq1_1", value = selectedSeq)
            fill = 1
        }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq2_1)==""){
                updateTextAreaInput(session, "piRNAseq2_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq3_1)==""){
                updateTextAreaInput(session, "piRNAseq3_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq4_1)==""){
                updateTextAreaInput(session, "piRNAseq4_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq5_1)==""){
                updateTextAreaInput(session, "piRNAseq5_1", value = selectedSeq)
                fill = 1
            }}
        
        if(fill == 0){
            if(as.character(input$piRNAseq6_1)==""){
                updateTextAreaInput(session, "piRNAseq6_1", value = selectedSeq)
                fill = 1
            }}
        
        ##Send custom message
        if(fill == 0){
            showNotification("Construct has already 6 sequences.")
        }
        }
        
        ##Second
        if(clust == "2"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_2)==""){
                    updateTextAreaInput(session, "piRNAseq1_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_2)==""){
                    updateTextAreaInput(session, "piRNAseq2_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_2)==""){
                    updateTextAreaInput(session, "piRNAseq3_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_2)==""){
                    updateTextAreaInput(session, "piRNAseq4_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_2)==""){
                    updateTextAreaInput(session, "piRNAseq5_2", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_2)==""){
                    updateTextAreaInput(session, "piRNAseq6_2", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 6 sequences.")
            }
        }
        
        ##Third
        if(clust == "3"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_3)==""){
                    updateTextAreaInput(session, "piRNAseq1_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_3)==""){
                    updateTextAreaInput(session, "piRNAseq2_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_3)==""){
                    updateTextAreaInput(session, "piRNAseq3_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_3)==""){
                    updateTextAreaInput(session, "piRNAseq4_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_3)==""){
                    updateTextAreaInput(session, "piRNAseq5_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_3)==""){
                    updateTextAreaInput(session, "piRNAseq6_3", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq7_3)==""){
                    updateTextAreaInput(session, "piRNAseq7_3", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 7 sequences.")
            }
        }
        
        
        ###Four
        if(clust == "4"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_4)==""){
                    updateTextAreaInput(session, "piRNAseq1_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_4)==""){
                    updateTextAreaInput(session, "piRNAseq2_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_4)==""){
                    updateTextAreaInput(session, "piRNAseq3_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_4)==""){
                    updateTextAreaInput(session, "piRNAseq4_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_4)==""){
                    updateTextAreaInput(session, "piRNAseq5_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_4)==""){
                    updateTextAreaInput(session, "piRNAseq6_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq7_4)==""){
                    updateTextAreaInput(session, "piRNAseq7_4", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq8_4)==""){
                    updateTextAreaInput(session, "piRNAseq8_4", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 8 sequences.")
            }
        }
        
        ###Five
        if(clust == "5"){
            if(fill == 0){
                if(as.character(input$piRNAseq1_5)==""){
                    updateTextAreaInput(session, "piRNAseq1_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq2_5)==""){
                    updateTextAreaInput(session, "piRNAseq2_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq3_5)==""){
                    updateTextAreaInput(session, "piRNAseq3_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq4_5)==""){
                    updateTextAreaInput(session, "piRNAseq4_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq5_5)==""){
                    updateTextAreaInput(session, "piRNAseq5_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq6_5)==""){
                    updateTextAreaInput(session, "piRNAseq6_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq7_5)==""){
                    updateTextAreaInput(session, "piRNAseq7_5", value = selectedSeq)
                    fill = 1
                }}
            
            if(fill == 0){
                if(as.character(input$piRNAseq8_5)==""){
                    updateTextAreaInput(session, "piRNAseq8_5", value = selectedSeq)
                    fill = 1
                }}
            
            ##Send custom message
            if(fill == 0){
                showNotification("Construct has already 8 sequences.")
            }
        }
    })
    
    ####Other functions###########
    ##Convert coordinates
    ConvCooTr2cDNA = function(posS, posE, ExonS, ExonE){
        if(length(posS) != length(posE)){return(c())}
        if((length(posS)<1) | (length(posE)<1)){return(c())}
        if((length(ExonS)<1) | (length(ExonE)<1)){return(cbind(posS,posE))}
        dists=ExonE - ExonS + 1
        if(sum(dists)==0){return(cbind(posS, posE))}
        resS=c()
        resE=c()
        for(i in 1:length(posS)){
            idx=which(ExonS < posS[i])
            if(length(idx) > 0){
                resS=append(resS,posS[i] - sum(dists[idx]))
            }else{
                resS=append(resS,posS[i])
            }
            idx=which(ExonS < posE[i])
            if(length(idx) > 0){
                resE=append(resE,posE[i] - sum(dists[idx]))
            }else{
                resE=append(resE,posE[i])
            }
        }
        return(cbind(resS, resE))
    }
    
    ##############################
    
    ###Create ApeFIle as pasteLines####
    PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,xtraComments,xtraLines, BinRevComp, organism){
        if (is.null(sequence)){return(NULL)}
        if(!is.character(sequence)){return(c())}
        if(length(patterns) < 1 ){return(c(paste(sequence)))}
        if(length(patterns) != length(FWDcolors)){return(c())}
        if(length(REVcolors) != length(FWDcolors)){return(c())}
        if(length(tooltips) != length(FWDcolors)){return(c())}
        if(length(BinRevComp)==0){BinRevComp=rep(FALSE,length(patterns))}
        
        ##Save Lines
        FileLines=c()
        FileLines=append(FileLines,paste("LOCUS",paste(locus_name,sep="",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
        FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
        FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
        FileLines=append(FileLines,paste("VERSION",".",sep="     "))
        FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
        FileLines=append(FileLines,paste("ORGANISM",organism,sep="     "))
        posipat=c()
        ##Match sequences
        for(i in 1:length(patterns)){
            stpos=start(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
            edpos=end(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
            if(length(stpos)>0){
                posipat=rbind(posipat, cbind(stpos,edpos,rep(tooltips[i],length(stpos)),rep(FWDcolors[i],length(stpos)),rep(REVcolors[i],length(stpos))))
            }
        }
        
        if(!(is.null(posipat))){
            colnames(posipat)=c("start","end","label","fwdc","revc")
        }
        
        if(xtraComments[1] != ""){
            FileLines=append(FileLines,paste("COMMENT",xtraComments,sep="     "))
        }
        
        if(!(is.null(posipat))){
            for(i in 1:length(patterns)){
                tempat=as.character(patterns[i])
                if(BinRevComp[i]){tempat=as.character(reverseComplement(DNAString(tempat)))}
                FileLines=append(FileLines,paste("COMMENT",paste(as.character(tooltips[i]),tempat),sep="     "))
            }
        }
        
        FileLines=append(FileLines,paste("COMMENT","Generated using wormbuilder.org/piRNAi/",sep="     "))
        FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
        
        if(!(is.null(posipat))){
            FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
            for(n in 1:nrow(posipat)){
                nixt= which(posipat[n,3] == tooltips)
                if(BinRevComp[nixt]){
                    xnoteA="complement("
                    xnoteB=")"
                    }else{
                        xnoteA=""
                        xnoteB=""
                    }
                FileLines=append(FileLines,paste("     primer_bind     ",xnoteA,c(posipat[n,1]),"..",c(posipat[n,2]),xnoteB,"",sep=""))
                FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
                FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
                FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"",c(posipat[n,4]),"\"",sep=""))
                FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"",c(posipat[n,5]),"\"",sep=""))
                FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
            }
            
        }
        
        if(xtraLines != ""){
            FileLines=append(FileLines,readLines(xtraLines))
        }
        
        FileLines=append(FileLines,paste("ORIGIN"))
        
        Compseq=unlist(strsplit(sequence,""))
        
        partseq=c()
        
        for(i in seq(1,length(Compseq),10)){
            endseq=i+9
            if(length(Compseq)-i < 9){endseq=length(Compseq)}
            partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
            
        }
        
        i=1
        for(num in seq(1,length(Compseq),60)){
            index=as.character(num)
            spaces=paste(rep(" ",6-nchar(index)),collapse="")
            endseq=i+5
            if((length(partseq)-i) < 5){endseq=length(partseq)}
            FileLines=append(FileLines , paste(spaces,index," ",paste(partseq[i:(endseq)],collapse=" "),sep=""))
            
            i=i+6
        }
        
        FileLines=append(FileLines,paste("//"))
        
        return(FileLines)
    }
    
    ##Retrieve output ape
    output$DownApeOut <- downloadHandler(
        filename <- function() {
            paste("piRNAi", "gb", sep=".")
        },
        
        content <- function(file) {
            file.copy(paste("WorkingSpace/users/",session_id,"/piRNAs.txt", sep=""), file)
        }
    )
    
    ##Retrieve construct ape
    output$DownConOut <- downloadHandler(
        filename <- function() {
            paste("piRNAi_construct", "gb", sep=".")
        },
        
        content <- function(file) {
            file.copy(paste("WorkingSpace/users/",session_id,"/construct.txt", sep=""), file)
        }
    )
    
    ##Calculate GC
    CalculateGC = function(x){
        if(!is.character(x)){return(c())}
        x=toupper(x)
        vecseq=unlist(strsplit(x,""))
        return((countPattern("C",x)+countPattern("G",x))/length(vecseq))
    }
    
    ##Clean boxes
    observeEvent(input$actionclean, {
        output$downloadconstruct <- renderUI({})
        output$AdvancedFragment <- renderText({})
        
        clust=as.character(input$clustercon)
        if(clust == "1"){
        updateTextAreaInput(session, "piRNAseq1_1", value = "")
        updateTextAreaInput(session, "piRNAseq2_1", value = "")
        updateTextAreaInput(session, "piRNAseq3_1", value = "")
        updateTextAreaInput(session, "piRNAseq4_1", value = "")
        updateTextAreaInput(session, "piRNAseq5_1", value = "")
        updateTextAreaInput(session, "piRNAseq6_1", value = "")
        }
        if(clust == "2"){
            updateTextAreaInput(session, "piRNAseq1_2", value = "")
            updateTextAreaInput(session, "piRNAseq2_2", value = "")
            updateTextAreaInput(session, "piRNAseq3_2", value = "")
            updateTextAreaInput(session, "piRNAseq4_2", value = "")
            updateTextAreaInput(session, "piRNAseq5_2", value = "")
            updateTextAreaInput(session, "piRNAseq6_2", value = "")
        }
        if(clust == "3"){
            updateTextAreaInput(session, "piRNAseq1_3", value = "")
            updateTextAreaInput(session, "piRNAseq2_3", value = "")
            updateTextAreaInput(session, "piRNAseq3_3", value = "")
            updateTextAreaInput(session, "piRNAseq4_3", value = "")
            updateTextAreaInput(session, "piRNAseq5_3", value = "")
            updateTextAreaInput(session, "piRNAseq6_3", value = "")
            updateTextAreaInput(session, "piRNAseq7_3", value = "")
        }
        if(clust == "4"){
            updateTextAreaInput(session, "piRNAseq1_4", value = "")
            updateTextAreaInput(session, "piRNAseq2_4", value = "")
            updateTextAreaInput(session, "piRNAseq3_4", value = "")
            updateTextAreaInput(session, "piRNAseq4_4", value = "")
            updateTextAreaInput(session, "piRNAseq5_4", value = "")
            updateTextAreaInput(session, "piRNAseq6_4", value = "")
            updateTextAreaInput(session, "piRNAseq7_4", value = "")
            updateTextAreaInput(session, "piRNAseq8_4", value = "")
        }
        if(clust == "5"){
            updateTextAreaInput(session, "piRNAseq1_5", value = "")
            updateTextAreaInput(session, "piRNAseq2_5", value = "")
            updateTextAreaInput(session, "piRNAseq3_5", value = "")
            updateTextAreaInput(session, "piRNAseq4_5", value = "")
            updateTextAreaInput(session, "piRNAseq5_5", value = "")
            updateTextAreaInput(session, "piRNAseq6_5", value = "")
            updateTextAreaInput(session, "piRNAseq7_5", value = "")
            updateTextAreaInput(session, "piRNAseq8_5", value = "")
        }
        })
})  
