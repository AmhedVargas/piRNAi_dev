########piRNA User Interface####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###UI

#Load libraries
library(shiny)
library(shinythemes)
library(DT)

##Load data to display for common fluorophores
fluo=read.table("WorkingSpace/Fluo_sizes.tab", sep="\t", stringsAsFactors=F)
colnames(fluo)=c("Name","File","Size")

# Define User interface
shinyUI(
    fluidPage(
        tags$head(
            tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
        ),
        ##Custom extra styles: single sliders background and title of navbar 
        tags$style(type = 'text/css', 
                   ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               "),
        tags$style(type='text/css', '#AdvancedFragment {white-space: pre-wrap;}'),
        tags$style(type='text/css', '#SimpleFragment {white-space: pre-wrap;}'),
        tags$style(type='text/css', '#piBoxes .form-group {margin-bottom: 0px; margin-top: 0px;}'),
        tags$style(type='text/css', "
                   table {
  font-family: arial, sans-serif;
  border-collapse: collapse;
  width: 100%;
}

td, th {
  border: 1px solid #dddddd;
  text-align: left;
  padding: 8px;
}

tr:nth-child(even) {
  background-color: #dddddd;
}
                   "),
        #Main tab pages
        navbarPage(
            title=actionLink("link_to_tabpanel_title", HTML("<b>piRNAi</b>")),
            windowTitle="piRNAi app",
            id = "panels",
            
            tabPanel("Designer",
                     mainPanel(
                         h2("Design piRNAi fragments"),
                         HTML("<h4>Target <i>C. elegans</i> and <i>C. briggsae</i> genes</h4>"),
                         tabsetPanel(
                            tabPanel("Simple",
                                     br(),
                         textAreaInput("geneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         actionButton("actiongenesearch", label = "Search"),
                         hr(),
                         uiOutput("DesignControls"),
                         hr(),
                         verbatimTextOutput("ErrorMessage"),
                         #tableOutput(otherPis),
                         htmlOutput("SelPiTabSummary"),
                         tableOutput('SelPiTab'),
                         verbatimTextOutput("SimpleFragment"),
                         uiOutput("downloadseq")
                         ),
                         tabPanel("Advanced",
                         h3("Select individual synthetic piRNAs"),
                    fluidRow(
                        ##Chose cluster
                        
                        radioButtons("clustercon", label = HTML("Select piRNA cluster
                                                     [<a href=\"\" onclick=\"$('#explain_cluster_adv').toggle(); return false;\">info</a>]
                                                     "),
                                     choices = list("Standard piRNA transgene (ClusterE, six piRNAs)" = 1, "ClusterG (six piRNAs)" = 2, "ClusterO (seven piRNAs)" = 3, "ClusterF (eigth piRNAs)" = 4, "New Cluster on 21ur-5764 (eigth piRNAs)" = 5), selected = 1, width='100%', inline= TRUE),
                        HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_cluster_adv\">
            We recommend to use the standard piRNAi clusterE (centered on the 21ur-1224 piRNA loci).
                                     </div></p>
                     "),
                        conditionalPanel(condition = "input.clustercon==1",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==2",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==3",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==4",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7"),
                        textAreaInput("piRNAseq8_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 8")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==5",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7"),
                        textAreaInput("piRNAseq8_5", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 8")),
                        br()
                        )),
                        column(8,
                        verbatimTextOutput("AdvancedFragment"),
                        uiOutput("downloadconstruct"))
                        ),
                        
                     fluidRow(
                         column(2,actionButton("actionclean", label = "Reset")),
                         column(2,actionButton("actionconstruct", label = "Generate piRNAi cluster"))),
                         verbatimTextOutput("AdvancedErrorMessage"),
                         hr(),
                         textAreaInput("Advancedgeneinput", label = "Gene target", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         actionButton("actionAdvsearch", label = "Search piRNAs"),
                        br(),
            br(),
            selectInput("FluoInput", label = HTML("<b>Common fluorophore
                    [<a href=\"\" onclick=\"$('#explain_n_fluo').toggle(); return false;\">info</a>]</b>"), choices=c(as.character(fluo$Name)), selected = 1),
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_n_fluo\">
            We have produced guide piRNAs that specifically targets common flurophore tags used in <i>C. elegans</i> labs. Please see the Download tab for more info.
                                     </div></p>
                     "),
            actionButton("ActionFluoSearch", label = "Show piRNAs"),
                        ##Dropdown menu
                        ##New button 
                         hr(),
                         uiOutput("AdvDesignControls"),
                         hr(),
                         DT::dataTableOutput('AllPiTab')
                     )
            ))),
            ###Background
            tabPanel("Background",
                     mainPanel(
                         h3("A database and cloning system for a genome-wide piRNA interference collection"),
                         HTML("<p align=\"justify\">
Recently, our lab has developed methods to silence genes via piRNAs instead (Priyadarshini <i>et al.</i>, 2021). piRNAs are a class of small RNAs that are active in the germline. We only recently learned how the piRNAs actually recognize genes (Heng-Chi lab paper and Mello lab paper). 
Knowing the rules, mean that we can design piRNAs to target specific genes and so, that is what we have done here.

This app helps to "),
                         actionLink("link_to_tabpanel_design", "design"),
                         HTML(" piRNAi fragments easily. Alternatively, you can "),
                         actionLink("link_to_tabpanel_downloads", "download"),
                         HTML(" our designs and see them in a genome browser."),
                     HTML(" </p>")
                     )
            ),
            ###About
            tabPanel("Downloads",
                     mainPanel(
                         h3("Tracks"),
                         HTML("<p align=\"justify\">
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Celegans.tar.gz\"><i>C. elegans</i> unique synthetic piRNAs (20-mers)</a><br>
                         <br><h3>Sequences</h3>
                         <table>
  <tr>
    <th>Name</th>
    <th>Lab</th>
    <th>GB file</th>
  </tr>
  <tr>
    <td>ce-GFP</td>
    <td>Froekjaer-Jensen lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/ce-GFP_AlJohani.ape\">ce-GFP_AlJohani</a></td>
  </tr>
    <tr>
    <td>ce-tagRFP</td>
    <td>Froekjaer-Jensen lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/ce-TagRFP-T_AlJohani.ape\">ce-TagRFP-T_AlJohani</a></td>
  </tr>
    <tr>
    <td>GFP</td>
    <td>Fire lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/gfp(S65C_Fire).ape\">gfp(S65C_Fire)</a></td>
  </tr>
    <tr>
    <td>mCherry</td>
    <td>Oegema lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mCherry_Oegema.ape\">mCherry_Oegema</a></td>
  </tr>
    <tr>
    <td>mCherry</td>
    <td>Seydoux lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mCherry_Seydoux.ape\">mCherry_Seydoux</a></td>
  </tr>
    <tr>
    <td>mKate2</td>
    <td>Dickinson lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mKate2_Dickinson.ape\">mKate2_Dickinson</a></td>
  </tr>
    <tr>
    <td>mNeonGreen</td>
    <td>Dickinson lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mNeonGreen_Dickinson.ape\">mNeonGreen_Dickinson</a></td>
  </tr>
</table>
                      </p>")
                     )
            ),
            ###About
            tabPanel("About",
                     mainPanel(
                         h3("The piRNA-based RNA interference app (\"piRNAi\")"),
                         HTML("<p align=\"justify\">
                      This website was generated via custom modified css/html code running in R via the shiny library.
                 <br>All the templates, libraries, and programs used to produce this site are under the MIT and GNU licenses.
                    <br>While the rationale behind the app is described in the background section, its implementation can be found at <a href=\"https://www.researchgate.net/profile/Amhed_Vargas_Velazquez\">Amhed Missael Vargas Velazquez's</a> <a href=\"https://github.com/AmhedVargas\">Github account</a></p>"),
                         h3("The Laboratory of Synthetic Genome Biology"),
                         HTML("<p align=\"justify\">
                 The Laboratory of Synthetic Genome Biology is located in building 2 - level 3 (Ibn Al-Haytham – Above Spine) at King Abdullah University of Science and Technology (KAUST).
                 <br><i>Contact info</i>:<br>Christian-Froekjaer Jensen, Ph.D. 
                 <br>Assistant Professor of Bioscience
                 <br>Laboratory of Synthetic Genome Biology
                 <br>Email: <a href=\"mailto:cfjensen@kaust.edu.sa\">cfjensen@kaust.edu.sa</a>
                      </p>")
                     )
            )
        ),
        hr(),
        HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
        HTML("<a href=\"mailto:cfjensen@kaust.edu.sa\">Contact us!</a>")
        
    )
)
