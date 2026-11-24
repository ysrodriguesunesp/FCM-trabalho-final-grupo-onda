# =============================================================================
# FCM – Trabalho final (Grupo Onda)
# Painel interativo: HIV/AIDS e sífilis nos municípios de SP, 2019–2023
#
# Estrutura da pasta (tudo isto é enviado ao shinyapps.io):
#   app/
#   ├── app.R
#   └── dados/
#       ├── AIDS_HIV.xlsx
#       ├── Sifilis.xlsx
#       ├── Densidade demografica.xlsx
#       ├── Indice de urbanizacao.xlsx
#       └── malha_municipios_SP.geojson
#
# Publicação:
#   rsconnect::deployApp("caminho/da/pasta/app", appFiles = c("app.R", "dados"),
#                        appName = "fcm-grupo-onda")
# =============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(sf)
library(leaflet)

sf_use_s2(FALSE)

# ---------------------------------------------------------------------------
# 0. Parâmetros
# ---------------------------------------------------------------------------
ESCALA        <- 1e5            # as taxas das planilhas são por 100 mil habitantes
ANOS_PANDEMIA <- 2020:2021
ANOS_BASE     <- c(2019, 2022)  # linha de base não pandêmica (2023 repete 2022)
COD_SMA       <- 3550209        # São Miguel Arcanjo
COD_CAP       <- 3550308        # São Paulo (capital)
COR_PANDEMIA  <- "#C0392B"
CORES_DOENCA  <- c("HIV/AIDS" = "#215B8F", "Sífilis" = "#B55F29")

# ---------------------------------------------------------------------------
# 1. Leitura dos dados e da malha municipal (executa uma vez, ao iniciar o app)
# ---------------------------------------------------------------------------
ler <- function(arquivo, variavel) {
  d <- read_excel(file.path("dados", arquivo))
  names(d) <- c("MUN", "CD_MUN", sub("2024", "2021", names(d)[-(1:2)]))
  d |>
    pivot_longer(-c(MUN, CD_MUN), names_to = "ano", values_to = variavel) |>
    mutate(ano = as.integer(ano), CD_MUN = as.numeric(CD_MUN))
}

malha <- st_read("dados/malha_municipios_SP.geojson", quiet = TRUE) |> st_make_valid()
malha$CD_MUN   <- as.numeric(malha$id)
malha$area_km2 <- as.numeric(st_area(st_transform(malha, 5880))) / 1e6
malha_web <- malha |>                                     # geometria simplificada para o mapa
  st_transform(5880) |>
  st_simplify(dTolerance = 250, preserveTopology = TRUE) |>
  st_transform(4326) |>
  select(CD_MUN)

base <- ler("AIDS_HIV.xlsx", "aids") |>
  left_join(ler("Sifilis.xlsx", "sifilis") |> select(-MUN), by = c("CD_MUN", "ano")) |>
  left_join(ler("Densidade demografica.xlsx", "densidade") |> select(-MUN), by = c("CD_MUN", "ano")) |>
  left_join(ler("Indice de urbanizacao.xlsx", "urbanizacao") |> select(-MUN), by = c("CD_MUN", "ano")) |>
  left_join(st_drop_geometry(malha)[, c("CD_MUN", "area_km2")], by = "CD_MUN") |>
  mutate(populacao = densidade * area_km2) |>
  rename(municipio = MUN)

# Correção da permutação capital × São Miguel Arcanjo (HIV/AIDS, 2019–2022)
for (a in 2019:2022) {
  i_s <- base$CD_MUN == COD_SMA & base$ano == a
  i_c <- base$CD_MUN == COD_CAP & base$ano == a
  taxa_s  <- base$aids[i_s]
  if (taxa_s < 100) taxa_s <- taxa_s * 1000               # desfaz a divisão por 1.000 feita na planilha
  casos_s <- taxa_s * base$populacao[i_s] / ESCALA
  casos_c <- base$aids[i_c] * base$populacao[i_c] / ESCALA
  base$aids[i_c] <- casos_s * ESCALA / base$populacao[i_c]
  base$aids[i_s] <- casos_c * ESCALA / base$populacao[i_s]
}

# Suavização bayesiana empírica global (Marshall, 1991)
bayes_empirico <- function(casos, pop) {
  m  <- sum(casos) / sum(pop)
  r  <- casos / pop
  s2 <- sum(pop * (r - m)^2) / sum(pop) - m / mean(pop)
  A  <- max(s2, 0)
  w  <- A / (A + m / pop)
  (w * r + (1 - w) * m) * ESCALA
}

dados <- base |>
  pivot_longer(c(aids, sifilis), names_to = "doenca", values_to = "taxa") |>
  mutate(doenca = recode(doenca, aids = "HIV/AIDS", sifilis = "Sífilis"),
         casos  = round(taxa * populacao / ESCALA)) |>
  group_by(doenca, ano) |>
  mutate(taxa_eb = bayes_empirico(casos, populacao)) |>
  ungroup()

anos <- sort(unique(dados$ano))
porte_2022 <- base |>
  filter(ano == 2022) |>
  transmute(CD_MUN, pop_2022 = populacao,
            porte = cut(populacao, c(0, 1e4, 5e4, 1e5, 5e5, Inf),
                        labels = c("< 10 mil", "10–50 mil", "50–100 mil", "100–500 mil", "> 500 mil")))

# Resumo por município (média 2019–2023 e variação na pandemia)
resumo_mun <- dados |>
  group_by(doenca, CD_MUN, municipio) |>
  summarise(media_eb  = mean(taxa_eb),
            base_eb   = mean(taxa_eb[ano %in% ANOS_BASE]),
            pand_eb   = mean(taxa_eb[ano %in% ANOS_PANDEMIA]),
            serie     = paste0(ano, ": ", formatC(taxa, format = "f", digits = 1, decimal.mark = ","),
                               collapse = " | "),
            .groups = "drop") |>
  mutate(variacao = pand_eb - base_eb) |>
  left_join(porte_2022, by = "CD_MUN") |>
  group_by(doenca) |>
  mutate(posicao = rank(-media_eb, ties.method = "min")) |>
  ungroup()

# Taxa estadual ponderada com IC 95% (Poisson exato)
estado <- dados |>
  group_by(doenca, ano) |>
  summarise(casos = sum(casos), pop = sum(populacao), .groups = "drop") |>
  mutate(taxa   = casos / pop * ESCALA,
         ic_inf = qchisq(0.025, 2 * casos) / 2 / pop * ESCALA,
         ic_sup = qchisq(0.975, 2 * (casos + 1)) / 2 / pop * ESCALA)

# Razão observados / esperados na pandemia, por porte
razao_oe <- function(d) {
  b <- d |> filter(ano %in% ANOS_BASE)
  p <- d |> filter(ano %in% ANOS_PANDEMIA)
  esperados  <- sum(b$casos) / sum(b$populacao) * sum(p$populacao)
  observados <- sum(p$casos)
  tibble(observados, esperados, OE = observados / esperados,
         ic_inf = qchisq(0.025, 2 * observados) / 2 / esperados,
         ic_sup = qchisq(0.975, 2 * (observados + 1)) / 2 / esperados)
}
oe <- bind_rows(
  dados |> group_by(doenca) |> group_modify(~ razao_oe(.x)) |> mutate(estrato = "Estado de SP"),
  dados |> left_join(porte_2022, by = "CD_MUN") |> group_by(doenca, porte) |>
    group_modify(~ razao_oe(.x)) |> ungroup() |> mutate(estrato = paste("Porte:", porte)) |> select(-porte)
) |> ungroup()
oe$estrato <- factor(oe$estrato, levels = rev(unique(oe$estrato)))

nomes_mun <- sort(unique(dados$municipio))
fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")

sombra_pandemia <- function() {
  list(annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = COR_PANDEMIA, alpha = 0.08),
       annotate("text", x = 2020.5, y = Inf, label = "Pandemia", vjust = 1.5, colour = COR_PANDEMIA, fontface = "bold"))
}

# ---------------------------------------------------------------------------
# 2. Interface
# ---------------------------------------------------------------------------
ui <- page_sidebar(
  title = "HIV/AIDS e sífilis nos municípios de São Paulo (2019–2023)",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2E7D5B"),

  sidebar = sidebar(
    width = 300,
    selectInput("doenca", "Agravo", choices = c("HIV/AIDS", "Sífilis")),
    selectInput("ano", "Ano", choices = anos, selected = 2022),
    selectizeInput("municipio", "Buscar município", choices = c("", nomes_mun),
                   options = list(placeholder = "Digite o nome...")),
    accordion(
      open = FALSE,
      accordion_panel(
        "Opções da aba Relação",
        selectInput("xvar", "Variável explicativa",
                    choices = c("Densidade demográfica" = "densidade",
                                "Índice de urbanização" = "urbanizacao")),
        numericInput("limite", "Excluir taxas acima de (por 100 mil hab.)", value = 300, min = 0, step = 50),
        checkboxInput("excluir_zero", "Excluir municípios com taxa zero", FALSE)
      )
    ),
    helpText("Taxas por 100 mil habitantes. Detalhes metodológicos na aba 'Sobre os dados'.")
  ),

  layout_columns(
    col_widths = c(4, 4, 4),
    value_box("Taxa estadual no ano (por 100 mil)", textOutput("vb_taxa"), theme = "primary"),
    value_box("Casos notificados no ano", textOutput("vb_casos"), theme = "secondary"),
    value_box("Pandemia: observados / esperados", textOutput("vb_oe"), textOutput("vb_oe_ic"), theme = "info")
  ),

  navset_card_tab(
    nav_panel(
      "Mapa",
      layout_columns(
        col_widths = c(8, 4),
        div(
          radioButtons("medida", NULL, inline = TRUE,
                       choices = c("Taxa bruta do ano" = "taxa",
                                   "Taxa suavizada do ano" = "taxa_eb",
                                   "Média suavizada 2019–2023" = "media_eb",
                                   "Variação na pandemia" = "variacao")),
          leafletOutput("mapa", height = 560)
        ),
        div(
          h5(textOutput("titulo_mun")),
          plotOutput("serie_mun", height = 300),
          tableOutput("tabela_mun")
        )
      )
    ),
    nav_panel("Pandemia",
              layout_columns(col_widths = c(6, 6),
                             plotOutput("estado", height = 430),
                             plotOutput("oe", height = 430))),
    nav_panel("Ranking", tableOutput("ranking")),
    nav_panel("Distribuição", plotOutput("histograma", height = 430)),
    nav_panel("Relação", plotOutput("dispersao", height = 430)),
    nav_panel(
      "Sobre os dados",
      markdown(paste(c(
        "### Fontes dos dados",
        "",
        "- **Casos de HIV/AIDS e de sífilis:** Sistema de Informação de Agravos de Notificação (SINAN), consultado no",
        "  **TabNet/DATASUS**, por município de residência, 2019 a 2023.",
        "- **População e densidade demográfica:** Instituto Brasileiro de Geografia e Estatística (**IBGE**), Censo Demográfico 2022.",
        "- **Limites municipais e área territorial:** malha municipal do **IBGE** (645 municípios do estado de São Paulo).",
        "- **Índice de urbanização:** razão entre a área urbanizada e a área total do município.",
        "",
        "### Indicadores",
        "",
        "- **Taxa bruta:** casos notificados ÷ população × 100.000 (casos por 100 mil habitantes).",
        "- **Taxa estadual ponderada:** soma dos casos ÷ soma da população dos municípios × 100.000, com intervalo de",
        "  confiança de 95% (Poisson exato).",
        "- **Taxa suavizada:** estimador bayesiano empírico global (Marshall, 1991), que aproxima as taxas de municípios",
        "  pequenos da taxa estadual na proporção da sua instabilidade.",
        "- **Média suavizada 2019–2023:** média das taxas suavizadas anuais de cada município.",
        "- **Variação na pandemia:** média das taxas suavizadas de 2020–2021 menos a média de 2019 e 2022, em casos por",
        "  100 mil habitantes (valores negativos indicam menos notificações na pandemia).",
        "- **Razão observados/esperados (O/E):** casos notificados em 2020–2021 ÷ casos esperados pela taxa de 2019 e 2022.",
        "  Valores abaixo de 1 indicam menos notificações que o esperado.",
        "",
        "### Notas metodológicas",
        "",
        "- Período pandêmico considerado: 2020–2021. As variações refletem a **notificação** dos agravos, que depende da",
        "  testagem e do acesso aos serviços, e não necessariamente a transmissão.",
        "- A área municipal foi calculada na projeção cônica equivalente de Albers (SIRGAS 2000, EPSG:5880).",
        "- As associações entre taxas e indicadores territoriais são ecológicas (nível municipal) e não se aplicam a",
        "  indivíduos.",
        "",
        "### Ferramentas",
        "",
        "R (shiny, bslib, leaflet, sf, ggplot2, dplyr), Python (pandas, geopandas, folium) e Git/GitHub.",
        "Mapa de fundo: Esri World Gray Canvas. Código e dados:",
        "[github.com/ysrodriguesunesp/FCM-trabalho-final-grupo-onda](https://github.com/ysrodriguesunesp/FCM-trabalho-final-grupo-onda).",
        "",
        "*Trabalho final da disciplina Ferramentas Computacionais de Modelagem — Grupo Onda.*"
      ), collapse = "\n"))
    )
  )
)

# ---------------------------------------------------------------------------
# 3. Servidor
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ano_sel <- reactive(as.integer(input$ano))
  do_ano  <- reactive(dados |> filter(doenca == input$doenca, ano == ano_sel()))

  # ---- indicadores ----
  output$vb_taxa  <- renderText(fmt(estado$taxa[estado$doenca == input$doenca & estado$ano == ano_sel()], 2))
  output$vb_casos <- renderText(fmt(estado$casos[estado$doenca == input$doenca & estado$ano == ano_sel()], 0))
  output$vb_oe    <- renderText({
    l <- oe |> filter(doenca == input$doenca, estrato == "Estado de SP")
    fmt(l$OE, 2)
  })
  output$vb_oe_ic <- renderText({
    l <- oe |> filter(doenca == input$doenca, estrato == "Estado de SP")
    paste0("IC 95% ", fmt(l$ic_inf, 2), "–", fmt(l$ic_sup, 2), " · 2020–2021 vs. 2019 e 2022")
  })

  # ---- mapa ----
  dados_mapa <- reactive({
    r <- resumo_mun |> filter(doenca == input$doenca)
    a <- do_ano() |> select(CD_MUN, taxa, taxa_eb, casos)
    malha_web |> inner_join(r, by = "CD_MUN") |> inner_join(a, by = "CD_MUN")
  })

  output$mapa <- renderLeaflet({
    leaflet(options = leafletOptions(minZoom = 6)) |>
      addProviderTiles(providers$Esri.WorldGrayCanvas) |>   # fundo sem exigência de chave de API
      setView(lng = -48.6, lat = -22.3, zoom = 7)
  })

  observe({
    d <- dados_mapa()
    v <- d[[input$medida]]
    if (input$medida == "variacao") {
      lim <- quantile(abs(v), 0.98, na.rm = TRUE)
      pal <- colorNumeric("RdBu", domain = c(-lim, lim), na.color = "#CCCCCC")
      cor <- pal(pmax(pmin(v, lim), -lim))
      titulo <- "Pandemia − linha de base<br>(vermelho = redução)"
      legenda <- function(m) addLegend(m, "bottomright", pal = pal, values = c(-lim, lim), title = titulo, opacity = 0.9)
    } else {
      cortes <- unique(c(0, quantile(v[v > 0], c(0.25, 0.5, 0.75, 0.9, 0.98), na.rm = TRUE), max(v, na.rm = TRUE)))
      pal <- colorBin(if (input$doenca == "HIV/AIDS") "Blues" else "Oranges", domain = v, bins = cortes, na.color = "#CCCCCC")
      cor <- pal(v)
      titulo <- "por 100 mil hab."
      legenda <- function(m) addLegend(m, "bottomright", pal = pal, values = v, title = titulo, opacity = 0.9,
                                       labFormat = labelFormat(digits = 1))
    }
    rotulos <- sprintf(
      "<b>%s</b><br>População (2022): %s<br>%s %s: taxa bruta %s | suavizada %s | casos %s<br>
       Média suavizada 2019–2023: %s (%sº de 645)<br>Variação na pandemia: %s<br><small>Taxas brutas: %s</small>",
      d$municipio, fmt(d$pop_2022, 0), input$doenca, input$ano, fmt(d$taxa), fmt(d$taxa_eb), fmt(d$casos, 0),
      fmt(d$media_eb), d$posicao, fmt(d$variacao), d$serie) |> lapply(HTML)

    leafletProxy("mapa", data = d) |>
      clearShapes() |> clearControls() |>
      addPolygons(layerId = ~CD_MUN, fillColor = cor, fillOpacity = 0.8, color = "#777777", weight = 0.4,
                  label = rotulos, highlightOptions = highlightOptions(weight = 2.5, color = "black", bringToFront = TRUE)) |>
      legenda()
    destacar(input$municipio, mover = FALSE)
  })

  destacar <- function(nome, mover = TRUE) {
    proxy <- leafletProxy("mapa") |> clearGroup("destaque")
    if (is.null(nome) || nome == "") return(invisible())
    cod <- dados$CD_MUN[dados$municipio == nome][1]
    g <- malha_web[malha_web$CD_MUN == cod, ]
    proxy |> addPolylines(data = g, color = "black", weight = 3.5, opacity = 1, group = "destaque")
    if (mover) {
      bb <- unname(st_bbox(g))
      proxy |> flyToBounds(bb[1], bb[2], bb[3], bb[4], options = list(maxZoom = 10))
    }
  }

  observeEvent(input$municipio, destacar(input$municipio), ignoreInit = TRUE)
  observeEvent(input$mapa_shape_click, {
    nome <- dados$municipio[dados$CD_MUN == input$mapa_shape_click$id][1]
    updateSelectizeInput(session, "municipio", selected = nome)
  })

  # ---- painel do município selecionado ----
  output$titulo_mun <- renderText({
    if (input$municipio == "") "Clique em um município no mapa ou use a busca" else input$municipio
  })

  output$serie_mun <- renderPlot({
    req(input$municipio != "")
    s <- dados |> filter(doenca == input$doenca, municipio == input$municipio)
    e <- estado |> filter(doenca == input$doenca)
    ggplot() +
      sombra_pandemia() +
      geom_line(data = e, aes(ano, taxa, linetype = "Estado (ponderada)"), colour = "black") +
      geom_line(data = s, aes(ano, taxa, linetype = "Taxa bruta"), colour = CORES_DOENCA[input$doenca], linewidth = 1) +
      geom_point(data = s, aes(ano, taxa), colour = CORES_DOENCA[input$doenca], size = 2.5) +
      geom_line(data = s, aes(ano, taxa_eb, linetype = "Taxa suavizada"), colour = CORES_DOENCA[input$doenca], alpha = 0.6) +
      scale_linetype_manual(values = c("Estado (ponderada)" = "dotted", "Taxa bruta" = "solid", "Taxa suavizada" = "dashed"),
                            name = NULL) +
      scale_y_continuous(limits = c(0, NA)) +
      labs(x = NULL, y = "Taxa por 100 mil hab.", title = input$doenca) +
      theme_minimal(base_size = 13) + theme(legend.position = "bottom")
  })

  output$tabela_mun <- renderTable({
    req(input$municipio != "")
    dados |>
      filter(doenca == input$doenca, municipio == input$municipio) |>
      transmute(Ano = as.character(ano), `População` = fmt(populacao, 0), Casos = fmt(casos, 0),
                `Taxa bruta` = fmt(taxa), `Taxa suavizada` = fmt(taxa_eb))
  })

  # ---- pandemia ----
  output$estado <- renderPlot({
    e <- estado |> filter(doenca == input$doenca)
    ggplot(e, aes(ano, taxa)) +
      sombra_pandemia() +
      geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), fill = CORES_DOENCA[input$doenca], alpha = 0.2) +
      geom_line(colour = CORES_DOENCA[input$doenca], linewidth = 1.2) +
      geom_point(colour = CORES_DOENCA[input$doenca], size = 3) +
      geom_text(aes(label = fmt(taxa)), vjust = -1.2, size = 3.8) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
      labs(x = NULL, y = "Taxa por 100 mil hab.",
           title = paste(input$doenca, "— taxa estadual ponderada (IC 95%)"),
           caption = "Variações refletem notificação, não necessariamente transmissão.") +
      theme_minimal(base_size = 13)
  })

  output$oe <- renderPlot({
    o <- oe |> filter(doenca == input$doenca)
    ggplot(o, aes(OE, estrato)) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_errorbarh(aes(xmin = ic_inf, xmax = ic_sup), height = 0.25, colour = CORES_DOENCA[input$doenca]) +
      geom_point(size = 3.5, colour = CORES_DOENCA[input$doenca]) +
      geom_text(aes(label = fmt(OE, 2)), vjust = -1, size = 3.6) +
      labs(x = "Observados / esperados em 2020–2021 (IC 95%)", y = NULL,
           title = "Pandemia vs. linha de base",
           caption = "< 1: menos notificações que o esperado; > 1: mais.") +
      theme_minimal(base_size = 13)
  })

  # ---- ranking ----
  output$ranking <- renderTable({
    resumo_mun |>
      filter(doenca == input$doenca) |>
      arrange(posicao) |>
      slice_head(n = 20) |>
      transmute(`Posição` = posicao, `Município` = municipio, `População (2022)` = fmt(pop_2022, 0),
                `Média suavizada 2019–2023` = fmt(media_eb), `Variação na pandemia` = fmt(variacao),
                `Taxas brutas por ano` = serie)
  })

  # ---- distribuição e relação (com filtros de sensibilidade) ----
  filtrado <- reactive({
    d <- do_ano() |> filter(taxa <= input$limite)
    if (input$excluir_zero) d <- d |> filter(taxa > 0)
    d
  })

  output$histograma <- renderPlot({
    ggplot(do_ano(), aes(taxa)) +
      geom_histogram(bins = 40, fill = CORES_DOENCA[input$doenca], colour = "white") +
      labs(x = paste(input$doenca, "– taxa por 100 mil habitantes"), y = "Número de municípios",
           title = paste("Distribuição das taxas municipais em", input$ano)) +
      theme_minimal(base_size = 14)
  })

  output$dispersao <- renderPlot({
    d <- filtrado()
    rho <- suppressWarnings(cor(d$taxa, d[[input$xvar]], method = "spearman"))
    g <- ggplot(d, aes(.data[[input$xvar]], taxa)) +
      geom_point(alpha = 0.5, colour = CORES_DOENCA[input$doenca]) +
      geom_smooth(method = "loess", formula = y ~ x, colour = COR_PANDEMIA) +
      labs(x = if (input$xvar == "densidade") "Densidade demográfica (hab/km², escala log)" else "Índice de urbanização",
           y = paste(input$doenca, "– taxa por 100 mil hab."),
           title = paste0("Taxa municipal × ", tolower(names(which(c("Densidade demográfica" = "densidade",
                                                                   "Índice de urbanização" = "urbanizacao") == input$xvar))),
                          " em ", input$ano, "  |  ρ de Spearman = ", fmt(rho, 2), "  (n = ", nrow(d), ")")) +
      theme_minimal(base_size = 14)
    if (input$xvar == "densidade") g <- g + scale_x_log10()
    g
  })
}

shinyApp(ui, server)
