" Is to be a local .vimrc ------------------------------------------------------

" The path for vim to search for files
set path+=../../../dealii-9.2.0/include

" English language
set spell spelllang=en_us

"Set fold method
set foldmethod=syntax

" To prevent Vim pretending it is Vi
set nocompatible

" Detecet file type of the file that is edited
filetype off

" Set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim

" Download plug-ins to the ~/.vim/plugged/ directory
call vundle#begin('~/.vim/plugged')

"------------------------------------------------------------------------------
" Let Vundle manage Vundle
Plugin 'VundleVim/Vundle.vim'

" Language packs. Language as C/C++, not like English.
Plugin 'sheerun/vim-polyglot'

" Source File Outliner
Plugin 'preservim/tagbar'

" YouCompleteMe
Plugin 'ycm-core/YouCompleteMe'

" Tags appear in the same order as in the file
let g:tagbar_sort = 0

" Tags are sorted by name
"let g:tagbar_sort = 1

"------------------------------------------------------------------------------

call vundle#end()

filetype plugin indent on

" Enable syntax highlighting
syntax enable

" Display all matching files when tab complete
set wildmenu

" Switches on line numbers
set number

" Wraps a line after 80 symbols
set textwidth=80

" Allow to use local .vimrc files
"set exrc

" Use mouse
set mouse=a

" Set a color scheme
colorscheme desert

" The following will make vim open new windows on the right side
set splitright

" Shows the line and column nubber in the status bar
set ruler

" The following commands set the spaces and tab characters visible
" set list
set listchars=tab:→\ ,space:·,nbsp:␣,trail:•,eol:¶,precedes:«,extends:»
set nolist

" Mapping function keys
map #5 :TagbarToggle <CR>
map #6 :YcmCompleter GoToDeclaration <CR>
map #7 :YcmCompleter GoToDefinition <CR>
map #8 :YcmCompleter GoToInclude <CR>
map #9 :YcmCompleter GetType <CR>

"Also very useful for jumping:
"Cttl-O, Ctrl-I see vim help

" Source File Outliner settings-------------------------------------------------

" Focus the panel when opening it
let g:tagbar_autofocus = 1

" Highlight the active tag
let g:tagbar_autoshowtag = 1

" Make panel vertical and place on the right
let g:tagbar_position = 'botright vertical'

" Tagbar's window width
let g:tagbar_width=85

" You complete me settings -----------------------------------------------------

" Mapping to close the completion menu (default <C-y>)
let g:ycm_key_list_stop_completion = ['<C-x>']

" Set filetypes where YCM will be turned on
let g:ycm_filetype_whitelist = { 'cpp':1, 'h':2, 'hpp':3, 'c':4, 'cxx':5 }

" Close preview window after completing the insertion
let g:ycm_autoclose_preview_window_after_insertion = 1
let g:ycm_autoclose_preview_window_after_completion = 1

" Don't confirm python conf
let g:ycm_confirm_extra_conf = 0

" Always populae diagnostics list
let g:ycm_always_populate_location_list = 1

" Enable line highligting diagnostics
let g:ycm_enable_diagnostic_signs = 1

" Open location list to view diagnostics
let g:ycm_open_loclist_on_ycm_diags = 1

" Max number of completion suggestions
let g:ycm_max_num_candidates = 20

" Max number of identifier-based suggestions
let g:ycm_max_num_identifier_candidates = 10

" Enable completion menu
let g:ycm_auto_trigger = 1

" Show diagnostic display features
let g:ycm_show_diagnostic_ui = 1

" The error symbol in Vim gutter
let g:ycm_error_symbol = '!'

" Display icons in Vim's gutter, error, warnings
let g:ycm_enable_diagnostic_signs = 1

" Highlight regions of diagnostic text
let g:ycm_enable_diagnostic_highlighting = 1

" Echo line's diagnostic that cursor is on
let g:ycm_echo_current_diagnostic = 1

"highlight Comment ctermfg=lightblue
"highlight Pmenu ctermfg=2 ctermbg=3 guifg=#ffffff guibg=#000000

highlight Pmenu ctermfg=White ctermbg=DarkGray guifg=#ffffff guibg=#000000

:hi SpellBad term=reverse cterm=bold ctermbg=9 ctermfg=0 gui=undercurl guisp=Red

