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

###Extra javascript function for input text
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

# Define User interface
shinyUI(
    fluidPage(
        tags$head(
            tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
        ),
        ##Javascript code
        tags$head(tags$script(HTML(jscode))),
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
                                     #tagAppendAttributes(
                         textAreaInput("geneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         #`data-proxy-click` = "actiongenesearch"),
                         actionButton("actiongenesearch", label = "Search"),
                         hr(),
                         uiOutput("DesignControls"),
                         hr(),
                         verbatimTextOutput("ErrorMessage"),
                         #tableOutput(otherPis),
                         htmlOutput("SelPiTabSummary"),
                         tableOutput('SelPiTab'),
                         verbatimTextOutput("SimpleFragment"),
                         htmlOutput("RefMesSim"),
                         uiOutput("downloadseq")
                         ),
                         tabPanel("Advanced",
                         h3("Select individual synthetic guide piRNAs"),
                    fluidRow(
                        ##Chose cluster
                        
                        radioButtons("clustercon", label = HTML("Select piRNA cluster
                                                     [<a href=\"\" onclick=\"$('#explain_cluster_adv').toggle(); return false;\">info</a>]
                                                     "),
                                     choiceNames = list(HTML("Standard piRNA transgene (Cluster<sub>E</sub>, six sg-piRNAs)"), HTML("Cluster<sub>G</sub> (six sg-piRNAs)"), HTML("Cluster<sub>O</sub> (seven sg-piRNAs)"), HTML("Cluster<sub>F</sub> (eight sg-piRNAs)")),choiceValues=c(1,2,3,4), selected = 1, width='100%', inline= TRUE),
                        HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_cluster_adv\">
            We recommend using the standard piRNAi Cluster<sub>E</sub> (centered on the <i>21ur-1224</i> piRNA loci). <b>Note:</b> Cluster<sub>O</sub> and Cluster<sub>F</sub> were tested with six sg-piRNAs but can accomodate additional sg-piRNAs, which may improve silencing.
                                     </div></p>
                     "),
                        conditionalPanel(condition = "input.clustercon==1",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_1", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 1"),
                        textAreaInput("piRNAseq2_1", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 2"),
                        textAreaInput("piRNAseq3_1", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 3"),
                        textAreaInput("piRNAseq4_1", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 4"),
                        textAreaInput("piRNAseq5_1", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 5"),
                        textAreaInput("piRNAseq6_1", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==2",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_2", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 1"),
                        textAreaInput("piRNAseq2_2", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 2"),
                        textAreaInput("piRNAseq3_2", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 3"),
                        textAreaInput("piRNAseq4_2", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 4"),
                        textAreaInput("piRNAseq5_2", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 5"),
                        textAreaInput("piRNAseq6_2", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==3",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 1"),
                        textAreaInput("piRNAseq2_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 2"),
                        textAreaInput("piRNAseq3_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 3"),
                        textAreaInput("piRNAseq4_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 4"),
                        textAreaInput("piRNAseq5_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 5"),
                        textAreaInput("piRNAseq6_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 6"),
                        textAreaInput("piRNAseq7_3", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 7")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==4",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 1"),
                        textAreaInput("piRNAseq2_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 2"),
                        textAreaInput("piRNAseq3_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 3"),
                        textAreaInput("piRNAseq4_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 4"),
                        textAreaInput("piRNAseq5_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 5"),
                        textAreaInput("piRNAseq6_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 6"),
                        textAreaInput("piRNAseq7_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 7"),
                        textAreaInput("piRNAseq8_4", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 8")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==5",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 1"),
                        textAreaInput("piRNAseq2_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 2"),
                        textAreaInput("piRNAseq3_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 3"),
                        textAreaInput("piRNAseq4_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 4"),
                        textAreaInput("piRNAseq5_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 5"),
                        textAreaInput("piRNAseq6_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 6"),
                        textAreaInput("piRNAseq7_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 7"),
                        textAreaInput("piRNAseq8_5", label="", rows=1, cols=21, resize="none", placeholder = "sg-piRNA 8")),
                        br()
                        )),
                        column(8,
                        verbatimTextOutput("AdvancedFragment"),
                        htmlOutput("RefMesAdv"),
                        uiOutput("downloadconstruct"))
                        ),
                        
                     fluidRow(
                         column(2,actionButton("actionclean", label = "Reset")),
                         column(2,actionButton("actionconstruct", label = "Generate piRNAi cluster"))),
                         verbatimTextOutput("AdvancedErrorMessage"),
                         hr(),
                          #tagAppendAttributes(
                         textAreaInput("Advancedgeneinput", label = "Gene target", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         #`data-proxy-click` = "actionAdvsearch"),
                         actionButton("actionAdvsearch", label = "Search piRNAs"),
                        br(),
            br(),
            selectInput("FluoInput", label = HTML("<b>Common fluorophore
                    [<a href=\"\" onclick=\"$('#explain_n_fluo').toggle(); return false;\">info</a>]</b>"), choices=c(as.character(fluo$Name)), selected = 1),
            HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_n_fluo\">
            We have designed sg-piRNAs that specifically target common fluorophores used in <i>C. elegans</i> labs. Please see the Download tab for primary DNA sequences. 
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
#             ###Background
#             tabPanel("Background",
#                      mainPanel(
#                          h3("A database and cloning system for a genome-wide piRNA interference collection"),
#                          HTML("<p align=\"justify\">
# Recently, our lab has developed methods to silence genes via piRNAs instead (Priyadarshini <i>et al.</i>, 2021). piRNAs are a class of small RNAs that are active in the germline. We only recently learned how the piRNAs actually recognize genes (Heng-Chi lab paper and Mello lab paper). 
# Knowing the rules, mean that we can design piRNAs to target specific genes and so, that is what we have done here.
# 
# This app helps to "),
#                          actionLink("link_to_tabpanel_design", "design"),
#                          HTML(" piRNAi fragments easily. Alternatively, you can "),
#                          actionLink("link_to_tabpanel_downloads", "download"),
#                          HTML(" our designs and see them in a genome browser."),
#                      HTML(" </p>")
#                      )
#             ),
            ###About
            tabPanel("Downloads",
                     mainPanel(
                         h3("Tracks"),
                         HTML("<br>Bed tracks with <i>C. elegans</i>  (WS270) guide-piRNAs for each specifity criteria.
                              
                              <p align=\"justify\">
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/piRNAi_Cel_WS270.tar.gz\">Download (307 MB)</a><br>
                         <br><h3>Sequences</h3>
                         <table>
  <tr>
    <th>Name</th>
    <th>Lab</th>
    <th>GenBank file (ApE annotations)</th>
  </tr>
  <tr>
    <td>ce-GFP</td>
    <td>Froekjaer-Jensen lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/ce-GFP_AlJohani.ape\">ce-gfp.gb</a></td>
  </tr>
    <tr>
    <td>ce-tagRFP</td>
    <td>Froekjaer-Jensen lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/ce-TagRFP-T_AlJohani.ape\">ce-tagRFP-T.gb</a></td>
  </tr>
    <tr>
    <td>GFP</td>
    <td>Fire lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/gfp(S65C_Fire).ape\">gfp(S65C).gb</a></td>
  </tr>
    <tr>
    <td>mCherry</td>
    <td>Oegema lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mCherry_Oegema.ape\">mCherry(Oegema).gb</a></td>
  </tr>
    <tr>
    <td>mCherry</td>
    <td>Seydoux lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mCherry_Seydoux.ape\">mCherry(Seydoux).gb</a></td>
  </tr>
    <tr>
    <td>mKate2</td>
    <td>Dickinson lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mKate2_Dickinson.ape\">mKate2.gb</a></td>
  </tr>
    <tr>
    <td>mNeonGreen</td>
    <td>Dickinson lab</td>
    <td><a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/piRNAi_dev/Apes/mNeonGreen_Dickinson.ape\">mNeonGreen.gb</a></td>
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
                    <br>While the rationale behind the app is described in the background section, its implementation can be found at <a href=\"https://www.researchgate.net/profile/Amhed_Vargas_Velazquez\">Amhed Missael Vargas Velazquez</a> | <a href=\"https://github.com/AmhedVargas\">Github account</a></p>"),
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
