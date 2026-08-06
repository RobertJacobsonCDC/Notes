import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'Math and Models',
  description: 'Notes on mathematical statistics and statistical modeling',
  cleanUrls: true,
  markdown: {
    math: true,
    image: {
      lazyLoading: true,
    },
    config(md) {
      md.core.ruler.after('block', 'github-math-compat', (state) => {
        for (const token of state.tokens) {
          if (token.type === 'fence' && token.info.trim() === 'math') {
            token.type = 'math_block'
            token.tag = 'math'
            continue
          }

          if (token.type === 'inline') {
            token.content = token.content.replace(
              /\$`([^`\n]+)`\$/g,
              (_, math: string) => `$${math}$`,
            )
          }
        }
      })
    },
  },
  themeConfig: {
    nav: [{ text: 'Probability Theory', link: '/' }],
    sidebar: [
      {
        text: 'Probability Theory',
        items: [
          { text: 'Introduction', link: '/' },
          { text: 'Probability Spaces', link: '/01_probability_space' },
          { text: 'Modeling Wait Time', link: '/02_wait_time' },
          { text: 'Poisson Processes', link: '/03_poisson_process' },
        ],
      },
    ],
    outline: {
      level: [2, 3],
    },
    search: {
      provider: 'local',
    },
  },
})
